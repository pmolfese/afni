// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Test adapter for the prepared real-library curves, before matrix sampling.
//! Arguments: duration TR peak|library. Prints time rows by twenty candidates.
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.len() != 3 {
        return Err("usage: hrf_probe duration TR peak|library".into());
    }
    let curves =
        afni_xmat_maker::sampled::glmsingle(args[0].parse()?, args[1].parse()?, args[2] == "peak")?;
    for row in 0..curves[0].values().len() {
        for (i, curve) in curves.iter().enumerate() {
            if i > 0 {
                print!(" ");
            }
            print!("{:.16e}", curve.values()[row]);
        }
        println!();
    }
    Ok(())
}
