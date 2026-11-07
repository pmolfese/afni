// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! EXPR syntax and scalar evaluation, independent of AFNI libraries.
//!
//! Language/domain conventions are adapted from AFNI parser.c/parser.inc and
//! parser_int.c (Medical College of Wisconsin, 1994–2000, CC BY 4.0).
//! Unlike the original bytecode interpreter, this module builds a small owned
//! expression tree. `Box<Node>` gives recursive children a fixed-size pointer;
//! Rust automatically frees the entire tree when its owner is dropped.
//! Parsing happens once; evaluation borrows the tree for each response sample.
//!
//! Random calls use an explicitly documented local 48-bit stream, not the
//! process-wide AFNI clock/PID seed. Cell permits updating just that small
//! state through a shared reference; it is not shared between threads.

use crate::{
    error::{Error, Result},
    expression_math as math,
};
use std::cell::Cell;

#[derive(Debug)]
pub(crate) struct Expression {
    root: Node,
    random: Cell<u64>,
}
#[derive(Debug)]
enum Node {
    Number(f64),
    Variable(u8),
    Negate(Box<Node>),
    Binary(u8, Box<Node>, Box<Node>),
    Call(String, Vec<Node>),
}

impl Expression {
    /// AFNI's NI_decode_string_list splits at whitespace or '~'. Shell quotes
    /// around the WHOLE basis argument are removed by the shell before we run;
    /// quotes embedded inside that argument do not protect expression spaces.
    pub(crate) fn parse_list(text: &str) -> Result<Vec<Self>> {
        let words: Vec<&str> = text
            .split(|ch: char| ch.is_ascii_whitespace() || ch == '~')
            .filter(|word| !word.is_empty())
            .collect();
        if words.is_empty() || words.len() > 1024 {
            return Err(Error::input("EXPR requires 1..1024 expressions"));
        }
        words
            .iter()
            .map(|word| {
                let cleaned: String = word
                    .chars()
                    .filter(|c| !c.is_whitespace())
                    .flat_map(char::to_uppercase)
                    .collect();
                if !cleaned.is_ascii() || cleaned.len() > 4096 {
                    return Err(Error::input(
                        "EXPR must be ASCII and at most 4096 bytes per component",
                    ));
                }
                let mut parser = Parser {
                    text: cleaned.as_bytes(),
                    at: 0,
                    depth: 0,
                    coordinate: false,
                };
                let root = parser.sum()?;
                if parser.at != parser.text.len() {
                    return Err(parser.error("unexpected trailing token"));
                }
                if cleaned != "1" && !parser.coordinate {
                    return Err(
                        parser.error("requires t, x, or z (only literal 1 may be constant)")
                    );
                }
                Ok(Self {
                    root,
                    random: Cell::new(0x1234_abcd_330e),
                })
            })
            .collect()
    }
    /// Reset this expression's private stream using srand48's 32-bit seed
    /// convention. Interior mutability is confined to this one integer.
    pub(crate) fn set_seed(&self, seed: u32) {
        self.random.set((u64::from(seed) << 16) | 0x330e);
    }
    pub(crate) fn evaluate(&self, t: f64, x: f64, z: f64) -> Option<f64> {
        self.root.evaluate([t, x, z], &self.random)
    }
}

/// Recursive descent mirrors the precedence levels: sum -> product -> unary
/// -> power -> atom. Each level consumes only its own operators, leaving the
/// next delimiter for its caller. Power associates to the right; unary minus
/// is outside power, so -t^2 means -(t^2), not (-t)^2.
struct Parser<'a> {
    text: &'a [u8],
    at: usize,
    depth: usize,
    coordinate: bool,
}
impl Parser<'_> {
    fn error(&self, message: &str) -> Error {
        Error::input(format!("EXPR byte {}: {message}", self.at + 1))
    }
    fn take(&mut self, c: u8) -> bool {
        if self.text.get(self.at) == Some(&c) {
            self.at += 1;
            true
        } else {
            false
        }
    }
    fn sum(&mut self) -> Result<Node> {
        self.depth += 1;
        if self.depth > 64 {
            return Err(self.error("nesting exceeds 64 levels"));
        }
        let mut node = self.product()?;
        loop {
            let op = if self.take(b'+') {
                b'+'
            } else if self.take(b'-') {
                b'-'
            } else {
                break;
            };
            node = Node::Binary(op, Box::new(node), Box::new(self.product()?));
        }
        self.depth -= 1;
        Ok(node)
    }
    fn product(&mut self) -> Result<Node> {
        let mut node = self.unary()?;
        loop {
            let op = if self.text.get(self.at..self.at + 2) == Some(b"**") {
                break;
            } else if self.take(b'*') {
                b'*'
            } else if self.take(b'/') {
                b'/'
            } else {
                break;
            };
            node = Node::Binary(op, Box::new(node), Box::new(self.unary()?));
        }
        Ok(node)
    }
    fn unary(&mut self) -> Result<Node> {
        if self.take(b'-') {
            self.depth += 1;
            if self.depth > 64 {
                return Err(self.error("too many unary operators"));
            }
            let child = self.unary()?;
            self.depth -= 1;
            Ok(Node::Negate(Box::new(child)))
        } else if self.take(b'+') {
            self.depth += 1;
            if self.depth > 64 {
                return Err(self.error("too many unary operators"));
            }
            let child = self.unary()?;
            self.depth -= 1;
            Ok(child)
        } else {
            self.power()
        }
    }
    fn power(&mut self) -> Result<Node> {
        let base = self.atom()?;
        let power = if self.take(b'^') {
            true
        } else if self.text.get(self.at..self.at + 2) == Some(b"**") {
            self.at += 2;
            true
        } else {
            false
        };
        if power {
            self.depth += 1;
            if self.depth > 64 {
                return Err(self.error("power nesting exceeds 64"));
            }
            let exponent = self.unary()?;
            self.depth -= 1;
            Ok(Node::Binary(b'^', Box::new(base), Box::new(exponent)))
        } else {
            Ok(base)
        }
    }
    fn atom(&mut self) -> Result<Node> {
        if self.take(b'(') {
            let node = self.sum()?;
            if !self.take(b')') {
                return Err(self.error("expected ')'"));
            }
            return Ok(node);
        }
        let start = self.at;
        if self
            .text
            .get(self.at)
            .is_some_and(|c| c.is_ascii_digit() || *c == b'.')
        {
            while self
                .text
                .get(self.at)
                .is_some_and(|c| c.is_ascii_digit() || *c == b'.')
            {
                self.at += 1;
            }
            if self
                .text
                .get(self.at)
                .is_some_and(|c| *c == b'E' || *c == b'D')
            {
                self.at += 1;
                if self
                    .text
                    .get(self.at)
                    .is_some_and(|c| *c == b'+' || *c == b'-')
                {
                    self.at += 1;
                }
                while self.text.get(self.at).is_some_and(u8::is_ascii_digit) {
                    self.at += 1;
                }
            }
            let text = std::str::from_utf8(&self.text[start..self.at])
                .unwrap()
                .replace('D', "E");
            let value: f64 = text.parse().map_err(|_| self.error("invalid number"))?;
            if !value.is_finite() {
                return Err(self.error("nonfinite constant"));
            }
            return Ok(Node::Number(value));
        }
        while self
            .text
            .get(self.at)
            .is_some_and(|c| c.is_ascii_alphanumeric() || *c == b'_')
        {
            self.at += 1;
        }
        if start == self.at {
            return Err(self.error("expected number, variable, or function"));
        }
        let name = std::str::from_utf8(&self.text[start..self.at]).unwrap();
        if !self.take(b'(') {
            if name == "PI" {
                return Ok(Node::Number(std::f64::consts::PI));
            }
            if name.len() == 1 && name.as_bytes()[0].is_ascii_uppercase() {
                self.coordinate |= matches!(name, "T" | "X" | "Z");
                return Ok(Node::Variable(name.as_bytes()[0]));
            }
            return Err(self.error("unknown constant or missing function parentheses"));
        }
        let Some(arity) = arity(name) else {
            return Err(self.error("unknown function"));
        };
        // These names are present in parser.inc but unavailable in the pinned
        // scalar reference: some invoke QQQERR, others have no evaluator arm.
        // Refuse them explicitly instead of returning an invented zero/identity.
        if matches!(
            name,
            "AI" | "DAI" | "BI" | "DBI" | "I0" | "I1" | "K0" | "K1" | "ABSEXTREME"
        ) {
            return Err(self.error("function is unimplemented in the pinned AFNI EXPR evaluator"));
        }
        let mut args = Vec::new();
        if !self.take(b')') {
            loop {
                if args.len() >= 128 {
                    return Err(self.error("function has more than 128 arguments"));
                }
                args.push(self.sum()?);
                if self.take(b')') {
                    break;
                }
                if !self.take(b',') {
                    return Err(self.error("expected ',' or ')'"));
                }
            }
        }
        if (arity >= 0 && args.len() != arity as usize) || args.is_empty() {
            return Err(self.error("wrong function argument count"));
        }
        let minimum = match name {
            "PAIRMAX" | "PAIRMIN" | "MOFN" | "ORSTAT" | "CHOOSE" | "AMONGST" => 2,
            "WITHIN" => 3,
            _ => 1,
        };
        if args.len() < minimum
            || (matches!(name, "PAIRMAX" | "PAIRMIN") && args.len() % 2 != 0)
            || (name == "WITHIN" && args.len() != 3)
        {
            return Err(self.error("invalid variable-argument function shape"));
        }
        Ok(Node::Call(name.to_owned(), args))
    }
}

impl Node {
    fn evaluate(&self, coordinates: [f64; 3], rng: &Cell<u64>) -> Option<f64> {
        let value = match self {
            Self::Number(v) => *v,
            Self::Variable(c) => match c {
                b'T' => coordinates[0],
                b'X' => coordinates[1],
                b'Z' => coordinates[2],
                _ => 0.0,
            },
            Self::Negate(child) => -child.evaluate(coordinates, rng)?,
            Self::Binary(op, left, right) => {
                let a = left.evaluate(coordinates, rng)?;
                let b = right.evaluate(coordinates, rng)?;
                match op {
                    b'+' => a + b,
                    b'-' => a - b,
                    b'*' => a * b,
                    b'/' => {
                        if b == 0.0 {
                            0.0
                        } else {
                            a / b
                        }
                    }
                    // AFNI leaves the base unchanged for an invalid real power,
                    // including 0^0 and a negative base with fractional exponent.
                    _ => {
                        if a > 0.0 || (a != 0.0 && b == b.trunc()) {
                            a.powf(b)
                        } else {
                            a
                        }
                    }
                }
            }
            Self::Call(name, args) => {
                let values: Option<Vec<_>> = args
                    .iter()
                    .map(|arg| arg.evaluate(coordinates, rng))
                    .collect();
                function(name, &values?, rng)?
            }
        };
        value.is_finite().then_some(value)
    }
}

/// Local drand48-compatible integer recurrence. The fixed initial state is a
/// deliberate reproducibility policy; AFNI instead seeds its global stream
/// from wall clock and PID. Matrix censoring never resamples this stream.
fn uniform(rng: &Cell<u64>) -> f64 {
    let state = rng.get().wrapping_mul(0x5deece66d).wrapping_add(11) & ((1 << 48) - 1);
    rng.set(state);
    state as f64 / (1_u64 << 48) as f64
}
/// Sort a temporary owned copy: evaluating a statistic must not reorder
/// the argument values used elsewhere in the expression tree.
fn median(a: &[f64]) -> f64 {
    let mut a = a.to_vec();
    a.sort_by(f64::total_cmp);
    let n = a.len();
    if n.is_multiple_of(2) {
        (a[n / 2 - 1] + a[n / 2]) * 0.5
    } else {
        a[n / 2]
    }
}
fn boolean(x: bool) -> f64 {
    if x { 1.0 } else { 0.0 }
}
/// Evaluate a parsed call. Arity was checked before this point, so indexed
/// arguments exist. None means undefined/nonfinite output, not a numeric zero.
/// Historical domain repairs are explicit match arms to make them reviewable.
fn function(name: &str, a: &[f64], rng: &Cell<u64>) -> Option<f64> {
    let x = a[0];
    let y = a.get(1).copied().unwrap_or(0.0);
    let value = match name {
        "SIN" => x.sin(),
        "COS" => x.cos(),
        "TAN" => x.tan(),
        "SIND" => (x * 0.01745329251994).sin(),
        "COSD" => (x * 0.01745329251994).cos(),
        "TAND" => (x * 0.01745329251994).tan(),
        "ASIN" => {
            if x.abs() <= 1.0 {
                x.asin()
            } else {
                x
            }
        }
        "ACOS" => {
            if x.abs() <= 1.0 {
                x.acos()
            } else {
                x
            }
        }
        "ATAN" => x.atan(),
        "ATAN2" => x.atan2(y),
        "SINH" => {
            if x.abs() < 87.5 {
                x.sinh()
            } else {
                x
            }
        }
        "COSH" => {
            if x.abs() < 87.5 {
                x.cosh()
            } else {
                x
            }
        }
        "TANH" => x.tanh(),
        "ASINH" => x.asinh(),
        "ACOSH" => {
            if x >= 1.0 {
                x.acosh()
            } else {
                x
            }
        }
        "ATANH" => {
            if x.abs() < 1.0 {
                x.atanh()
            } else {
                x
            }
        }
        "EXP" => x.min(87.5).exp(),
        "LOG" => {
            if x != 0.0 {
                x.abs().ln()
            } else {
                0.0
            }
        }
        "LOG10" => {
            if x != 0.0 {
                x.abs().log10()
            } else {
                0.0
            }
        }
        "ABS" => x.abs(),
        "INT" => x.trunc(),
        "SQRT" => x.abs().sqrt(),
        "CBRT" => x.cbrt(),
        "MAX" => x.max(y),
        "MIN" => x.min(y),
        "MOD" => {
            if y != 0.0 {
                x - y * (x / y).trunc()
            } else {
                0.0
            }
        }
        "LOGCOSH" => {
            if x.abs() < 87.5 {
                x.cosh().ln()
            } else {
                x
            }
        }
        "ERF" => math::erf(x),
        "ERFC" => math::erfc(x),
        "GAMMA" => math::gamma(x),
        "QG" => math::normal_q(x),
        "QGINV" => math::normal_q_inverse(x),
        "J0" => math::bessel(x, 0, false),
        "J1" => math::bessel(x, 1, false),
        "Y0" => math::bessel(x, 0, true),
        "Y1" => math::bessel(x, 1, true),
        "BELL2" => {
            if x.abs() <= 0.5 {
                1.0 - x.abs() * 1.3333333333333333 * x.abs()
            } else if x.abs() <= 1.5 {
                (1.5 - x.abs()).powi(2) * 0.666666666666667
            } else {
                0.0
            }
        }
        "RECT" => boolean(x.abs() <= 0.5),
        "TENT" => (1.0 - x.abs()).max(0.0),
        "STEP" | "ISPOSITIVE" => boolean(x > 0.0),
        "ISNEGATIVE" => boolean(x < 0.0),
        "POSVAL" => x.max(0.0),
        "BOOL" | "NOTZERO" => boolean(x != 0.0),
        "NOT" | "ISZERO" => boolean(x == 0.0),
        "EQUALS" => boolean(x == y),
        "ASTEP" => boolean(x.abs() > y),
        "AND" => boolean(a.iter().all(|&v| v != 0.0)),
        "OR" => boolean(a.iter().any(|&v| v != 0.0)),
        "MOFN" => boolean(a[1..].iter().filter(|&&v| v != 0.0).count() as i64 >= x as i64),
        "AMONGST" => boolean(a[1..].contains(&x)),
        "WITHIN" => boolean(x >= y && x <= a[2]),
        "ARGNUM" => a.iter().filter(|&&v| v != 0.0).count() as f64,
        "ARGMAX" => {
            if a.iter().all(|&v| v == 0.0) {
                0.0
            } else {
                let mut best = 0;
                for i in 1..a.len() {
                    if a[i] > a[best] {
                        best = i;
                    }
                }
                (best + 1) as f64
            }
        }
        "PAIRMAX" | "PAIRMIN" => {
            let half = a.len() / 2;
            let mut best = 0;
            for i in 1..half {
                if (name == "PAIRMAX" && a[i] > a[best]) || (name == "PAIRMIN" && a[i] < a[best]) {
                    best = i;
                }
            }
            a[half + best]
        }
        "MEDIAN" => median(a),
        "MAD" => {
            let mid = median(a);
            median(&a.iter().map(|v| (v - mid).abs()).collect::<Vec<_>>())
        }
        "MEAN" => a.iter().sum::<f64>() / a.len() as f64,
        "STDEV" | "SEM" => {
            if a.len() < 2 {
                0.0
            } else {
                let mean = a.iter().sum::<f64>() / a.len() as f64;
                let var = a.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (a.len() - 1) as f64;
                (var / if name == "SEM" { a.len() as f64 } else { 1.0 }).sqrt()
            }
        }
        "ORSTAT" => {
            let mut values = a[1..].to_vec();
            values.sort_by(f64::total_cmp);
            let index = (x as i64).clamp(1, values.len() as i64);
            values[index as usize - 1]
        }
        "CHOOSE" => {
            let index = x as i64;
            if index < 1 || index >= a.len() as i64 {
                0.0
            } else {
                a[index as usize]
            }
        }
        "IFELSE" => {
            if x != 0.0 {
                y
            } else {
                a[2]
            }
        }
        "HMODE" | "LMODE" => {
            let mut values = a.to_vec();
            values.sort_by(f64::total_cmp);
            let mut best = values[0];
            let mut best_count = 0;
            let mut i = 0;
            while i < values.len() {
                let mut j = i + 1;
                while j < values.len() && values[j] == values[i] {
                    j += 1;
                }
                if j - i > best_count || (name == "HMODE" && j - i == best_count) {
                    best = values[i];
                    best_count = j - i;
                }
                i = j;
            }
            best
        }
        "MINABOVE" => a[1..]
            .iter()
            .copied()
            .filter(|&v| v > x)
            .reduce(f64::min)
            .unwrap_or(x),
        "MAXBELOW" => a[1..]
            .iter()
            .copied()
            .filter(|&v| v < x)
            .reduce(f64::max)
            .unwrap_or(x),
        "EXTREME" | "ABSEXTREME" => {
            // EXTREME compares each absolute candidate against a SIGNED
            // accumulator in parser.c. Preserve this historical behavior;
            // e.g. extreme(-8,-3,2) is 2, not -8. ABSEXTREME is rejected
            // by parsing because its reference opcode is truncated.
            let mut best = 0.0;
            for &v in a {
                if v.abs() > best {
                    best = v;
                }
            }
            if best == 0.0 { x } else { best }
        }
        "ZTONE" => {
            if x <= 0.0 {
                0.0
            } else if x >= 1.0 {
                1.0
            } else {
                0.50212657
                    * ((std::f64::consts::FRAC_PI_2 * (1.6 * x - 0.8)).tan().tanh() + 0.99576486)
            }
        }
        "PLEG" => math::legendre(x as i32, y),
        "HRFBK4" | "HRFBK5" => math::hrf(x, y, if name == "HRFBK4" { 4 } else { 5 }),
        "RHDDC2" => math::rhdd(a),
        "ACFWXM" => math::acf_width(a),
        "GAMP" | "GAMQ" => {
            let (p, q) = crate::basis_extended::peak_width_double(x, y)?;
            if name == "GAMP" { p } else { q }
        }
        "ISPRIME" => {
            if x <= 0.0 || x > 2_147_483_647.0 || (x - x.round_ties_even()).abs() > 0.001 {
                -1.0
            } else if x.round_ties_even() < 2.0 {
                0.0
            } else {
                let n = x.round_ties_even() as u64;
                let mut d = 2;
                let mut prime = true;
                while d * d <= n {
                    if n.is_multiple_of(d) {
                        prime = false;
                        break;
                    }
                    d += 1;
                }
                boolean(prime)
            }
        }
        "URAN" => x * uniform(rng),
        "IRAN" => ((x + 1.0) * uniform(rng)).trunc(),
        "ERAN" => -x * uniform(rng).max(f64::MIN_POSITIVE).ln(),
        "LRAN" => x * (1.0 / uniform(rng).max(f64::MIN_POSITIVE) - 1.0).ln(),
        // Box–Muller targets a Gaussian directly. AFNI mixes a twelve-uniform
        // approximation with its cached Box–Muller sampler; we intentionally
        // do not reproduce that global state or its clock/PID seed.
        "GRAN" => {
            x + y
                * (-2.0 * uniform(rng).max(f64::MIN_POSITIVE).ln()).sqrt()
                * (std::f64::consts::TAU * uniform(rng)).cos()
        }
        _ => return math::statistical(name, a),
    };
    value.is_finite().then_some(value)
}

/// Frozen dispatch arities from parser.inc. Keeping the table explicit makes
/// supported spellings reviewable and lets tests audit every source entry.
fn arity(name: &str) -> Option<i32> {
    match name {
        "SIN" => Some(1),
        "COS" => Some(1),
        "TAN" => Some(1),
        "ASIN" => Some(1),
        "ACOS" => Some(1),
        "ATAN" => Some(1),
        "ATAN2" => Some(2),
        "SINH" => Some(1),
        "COSH" => Some(1),
        "TANH" => Some(1),
        "ASINH" => Some(1),
        "ACOSH" => Some(1),
        "ATANH" => Some(1),
        "EXP" => Some(1),
        "LOG" => Some(1),
        "LOG10" => Some(1),
        "ABS" => Some(1),
        "INT" => Some(1),
        "SQRT" => Some(1),
        "MAX" => Some(2),
        "MIN" => Some(2),
        "AI" => Some(1),
        "DAI" => Some(1),
        "I0" => Some(1),
        "I1" => Some(1),
        "J0" => Some(1),
        "J1" => Some(1),
        "K0" => Some(1),
        "K1" => Some(1),
        "Y0" => Some(1),
        "Y1" => Some(1),
        "BI" => Some(1),
        "DBI" => Some(1),
        "ERF" => Some(1),
        "ERFC" => Some(1),
        "GAMMA" => Some(1),
        "QG" => Some(1),
        "QGINV" => Some(1),
        "BELL2" => Some(1),
        "RECT" => Some(1),
        "STEP" => Some(1),
        "BOOL" => Some(1),
        "AND" => Some(-1),
        "OR" => Some(-1),
        "MOFN" => Some(-1),
        "ASTEP" => Some(2),
        "SIND" => Some(1),
        "COSD" => Some(1),
        "TAND" => Some(1),
        "MEDIAN" => Some(-1),
        "FICO_T2P" => Some(4),
        "FICO_P2T" => Some(4),
        "FICO_T2Z" => Some(4),
        "FITT_T2P" => Some(2),
        "FITT_P2T" => Some(2),
        "FITT_T2Z" => Some(2),
        "FIFT_T2P" => Some(3),
        "FIFT_P2T" => Some(3),
        "FIFT_T2Z" => Some(3),
        "FIZT_T2P" => Some(1),
        "FIZT_P2T" => Some(1),
        "FIZT_T2Z" => Some(1),
        "FICT_T2P" => Some(2),
        "FICT_P2T" => Some(2),
        "FICT_T2Z" => Some(2),
        "FIBT_T2P" => Some(3),
        "FIBT_P2T" => Some(3),
        "FIBT_T2Z" => Some(3),
        "FIBN_T2P" => Some(3),
        "FIBN_P2T" => Some(3),
        "FIBN_T2Z" => Some(3),
        "FIGT_T2P" => Some(3),
        "FIGT_P2T" => Some(3),
        "FIGT_T2Z" => Some(3),
        "FIPT_T2P" => Some(2),
        "FIPT_P2T" => Some(2),
        "FIPT_T2Z" => Some(2),
        "ZTONE" => Some(1),
        "LMODE" => Some(-1),
        "HMODE" => Some(-1),
        "GRAN" => Some(2),
        "URAN" => Some(1),
        "IRAN" => Some(1),
        "ERAN" => Some(1),
        "LRAN" => Some(1),
        "ORSTAT" => Some(-1),
        "TENT" => Some(1),
        "MAD" => Some(-1),
        "ARGMAX" => Some(-1),
        "ARGNUM" => Some(-1),
        "NOTZERO" => Some(1),
        "ISZERO" => Some(1),
        "EQUALS" => Some(2),
        "ISPOSITIVE" => Some(1),
        "ISNEGATIVE" => Some(1),
        "MEAN" => Some(-1),
        "STDEV" => Some(-1),
        "SEM" => Some(-1),
        "PLEG" => Some(2),
        "CDF2STAT" => Some(5),
        "STAT2CDF" => Some(5),
        "PAIRMAX" => Some(-1),
        "PAIRMIN" => Some(-1),
        "AMONGST" => Some(-1),
        "CBRT" => Some(1),
        "RHDDC2" => Some(3),
        "HRFBK4" => Some(2),
        "HRFBK5" => Some(2),
        "POSVAL" => Some(1),
        "NOT" => Some(1),
        "MOD" => Some(2),
        "WITHIN" => Some(-1),
        "MINABOVE" => Some(-1),
        "MAXBELOW" => Some(-1),
        "EXTREME" => Some(-1),
        "ABSEXTREME" => Some(-1),
        "CHOOSE" => Some(-1),
        "IFELSE" => Some(3),
        "LOGCOSH" => Some(1),
        "ACFWXM" => Some(4),
        "GAMP" => Some(2),
        "GAMQ" => Some(2),
        "ISPRIME" => Some(1),
        _ => None,
    }
}
