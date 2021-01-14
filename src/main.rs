mod wfa;
mod tmp;

fn main() {
    // tmp::tmp();
    wfa_main();
}

fn wfa_main() {
    let mut line = String::new(); 
    match std::io::stdin().read_line(&mut line) {
        Ok(_n) => {
            // split string into two
            let mut iter = line.split_whitespace();
            let str1 = iter.next().unwrap_or("");
            let str2 = iter.next().unwrap_or("");
            if str1 == "" || str2 == "" {
                std::process::exit(1);
            } 
            let opts = wfa::WFAOpts {
                a: 0,
                x: 4,
                e: 2,
                o: 6,
            };
            let score = wfa::wfa::<i8>(&opts, str1.as_bytes(), str2.as_bytes());
            println!("{}", score);
        }
        Err(error) => {
            println!("error: {}", error);
        }
    }
}
// 
//    g a t a c a
//  g 0 1 2 3 4 5
//  A 1 0 1 2 3 4
//  G 2 1 0 1 2 3 
//  A 3 2 1 0 1 2 
//  T 4 3 2 1 0 1
//  A 5 4 3 2 1 0