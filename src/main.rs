mod wfa;
mod tmp;

fn main() {
    let opts = wfa::WFAOpts {
        a : 0,
        x : 4,
        o : 6,
        e : 2,
    };
    let t = "GATACA";
    let q = "GAGATA";
    // let t = "GATACA";
    // let q = "GATACA";
    wfa::wfa::<i8>(opts, t.as_bytes(), q.as_bytes());
}
// 
//    G A T A C A
//  G 0 1 2 3 4 5
//  A 1 0 1 2 3 4
//  G 2 1 0 1 2 3 
//  A 3 2 1 0 1 2 
//  T 4 3 2 1 0 1
//  A 5 4 3 2 1 0