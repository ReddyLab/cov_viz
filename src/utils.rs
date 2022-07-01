use std::time::Instant;

pub fn print_pct_complete(total: i64) -> impl FnMut(i64) -> () {
    let five_pct_rc = total as f64 * 0.05;
    let mut pct_rc = five_pct_rc;
    let start_time = Instant::now();

    move |current| {
        if (current as f64) > pct_rc {
            println!(
                "Percent complete: {:.0}%, ({:.0}s)",
                ((current as f64) / (total as f64)) * 100.0,
                start_time.elapsed().as_secs()
            );
            pct_rc = (current as f64) + five_pct_rc;
        }
    }
}
