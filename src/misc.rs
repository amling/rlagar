use chrono::Local;

pub fn debug_log(msg: impl AsRef<str>) {
    let msg = msg.as_ref();
    eprintln!("{} - {}", Local::now().format("%Y%m%d %H:%M:%S"), msg);
}

pub fn debug_time<T>(label: impl AsRef<str>, cb: impl FnOnce() -> T) -> T {
    let label = label.as_ref();
    let t0 = std::time::Instant::now();
    debug_log(format!("Starting {}...", label));
    let ret = cb();
    debug_log(format!("Finished {}: {:?}", label, t0.elapsed()));
    return ret;
}
