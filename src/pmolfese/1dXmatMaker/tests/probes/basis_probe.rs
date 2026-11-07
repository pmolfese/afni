// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Test adapter for individual waveform samples, including dm families before
//! their M4 timing-file integration. This example is not an extra public CLI.
//! It uses the same public waveform API available to future Rust callers.
use afni_xmat_maker::basis::Basis;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Collect owned argument strings so the parsed formula can borrow one for
    // the duration of this call. `?` returns a parse error to main's caller.
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 5 || args.len() > 6 {
        return Err("usage: basis_probe FORMULA START STEP COUNT [EVENT_DURATION]".into());
    }
    let basis = Basis::parse(&args[1])?;
    let start: f64 = args[2].parse()?;
    let step: f64 = args[3].parse()?;
    let count: usize = args[4].parse()?;
    if count > 1_000_000 {
        return Err("probe count is too large".into());
    }
    let duration: Option<f64> = args.get(5).map(|s| s.parse()).transpose()?;
    for row in 0..count {
        let time = start + row as f64 * step;
        for component in 0..basis.components() {
            let value = if let Some(duration) = duration {
                basis.sample_with_duration(time, duration)
            } else {
                basis.sample(component, time)
            }
            .ok_or("invalid or nonfinite waveform sample")?;
            print!("{}{value:.17}", if component == 0 { "" } else { " " });
        }
        println!();
    }
    Ok(())
}
