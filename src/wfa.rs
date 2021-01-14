use num;
use std::option::Option;
use std::collections::VecDeque;

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
#[cfg(target_arch = "x86")]
use std::arch::x86::*;


pub struct WFAOpts {
    pub a: u8,
    pub x: u8,
    pub o: u8,
    pub e: u8,
}

impl Default for WFAOpts {
    fn default() -> Self {
        return WFAOpts {
            a: 0,
            x: 1,
            o: 1,
            e: 1,
        }
    }
}

fn init_wf<T>(l: usize) -> Vec<T>
where
    T:  std::clone::Clone + num::FromPrimitive,
{
    return vec![T::from_i32(-1).unwrap(); l];
}

// T: integer. use larger ints for larger sequences
// q.len() and t.len() must each be smaller than 2^sizeof(T)
// https://doi.org/10.1093/bioinformatics/btaa777
pub fn wfa<T>(opts: &WFAOpts, t: &[u8], q: &[u8]) -> usize
    where 
        T: num::PrimInt + num::Signed + num::Zero + num::One + num::FromPrimitive + num::traits::Saturating + std::clone::Clone + std::fmt::Debug,
{
    let max_len = num::pow(2, std::mem::size_of::<T>()*8-1)-1;
    if q.len() > max_len || t.len() > max_len {
        panic!("query length {} or template length {} greater than max: {}. Please consider using a longer integer type", q.len(), t.len(), max_len);
    }
    let k_max : usize = t.len() + q.len() - 2;
    let k_zero : usize = q.len() - 1;
    let a_k : usize = match t.len() >= q.len() { // index of bottom right corner of DP matrix
        true => t.len() - q.len() + k_zero,
        false => q.len() - t.len() + k_zero,
    };
    let a_off : usize = t.len(); // furthest offset of longest diagonal
    // current wavefronts
    let mut wfs_s : [Vec<T>; 3] = [init_wf(k_max+1), init_wf(k_max+1), init_wf(k_max+1)]; // m, i, d at current score
    wfs_s[0][k_zero] = T::zero();
    // these are buffers of past WFs needed for each of the operations
    let mut wfs_m : VecDeque<Vec<T>> = VecDeque::with_capacity(opts.x as usize);
    let mut wfs_o : VecDeque<Vec<T>> = VecDeque::with_capacity((opts.o + opts.e) as usize);
    let mut wfs_ie : VecDeque<Vec<T>> = VecDeque::with_capacity(opts.e as usize);
    let mut wfs_de : VecDeque<Vec<T>> = VecDeque::with_capacity(opts.e as usize);
    let mut bounds = [[k_zero, k_zero], [k_zero+1, k_zero+1], [k_zero-1, k_zero-1]];
    let mut s : usize = 0; // score (must be positive)
    loop {
        // wf_extend(m[s], q, t);
        // extend wavefront by searching down diagonals for matches
        for k in bounds[0][0]..std::cmp::min(bounds[0][1], k_max) + 1 {
            if k > k_max {
                break;
            }
            let mut v : usize = (wfs_s[0][k].to_i32().unwrap() + (k_zero as i32 - k as i32)) as usize;
            let mut h : usize = wfs_s[0][k].to_usize().unwrap_or(0);
            while v < q.len() && h < q.len() && q[v] == t[h] {
                wfs_s[0][k] = wfs_s[0][k] + T::one();
                v = v + 1;
                h = h + 1;
            }
        }
        // check to see if wavefront extended all the way to bottom right of DP matrix
        if wfs_s[0][a_k] >= T::from_usize(a_off).unwrap_or(T::zero()) {
            break;
        }
        // accumulate wfs in buffers (one for each kind of score)
        wfs_m.push_back(wfs_s[0].clone());
        wfs_o.push_back(wfs_s[0].clone());
        wfs_ie.push_back(wfs_s[1].clone());
        wfs_de.push_back(wfs_s[2].clone());


        // only increment by opts.e if gap open already happened
        s += match s >= (opts.o + opts.e + opts.e) as usize {
            true =>  *[opts.o + opts.e, opts.x, opts.e].iter().min().unwrap() as usize,
            false => *[opts.o + opts.e, opts.x].iter().min().unwrap() as usize,
        };

        // wf_next(m, i, d, q, t, s);
        // get bounds for new wavefront
        // wavefronts for [mismatch, ins, del]
        // set bounds of M
        let mut wf_m : Option<Vec<T>> = None;
        let mut wf_o : Option<Vec<T>> = None;
        let mut wf_ie : Option<Vec<T>> = None;
        let mut wf_de : Option<Vec<T>> = None;
        if s >= opts.x as usize {
            if let Some(i) = wfs_m.pop_front() {
                wf_m = Some(i);
                bounds[0][0] = bounds[0][0] - 1;
                bounds[0][1] = bounds[0][1] + 1;
            };
        }
        if s >= (opts.o + opts.e) as usize {
            wf_o =  wfs_o.pop_front() ;
        }
        // extend bounds of I and D only after one open+ext is done
        if s >= (opts.o + opts.e + opts.e) as usize { // when extension is possible
            if let Some(i) = wfs_ie.pop_front() {
                wf_ie = Some(i); 
                bounds[1][0] = bounds[1][0] - 1;
                bounds[1][1] = bounds[1][1] + 1;
            };
            if let Some(i) = wfs_de.pop_front() {
                wf_de = Some(i);
                bounds[2][0] = bounds[2][0] - 1;
                bounds[2][1] = bounds[2][1] + 1;
            };
        }
        // update I, D, M here.
        for k in bounds[1][0]..bounds[1][1]+1 {
            wfs_s[1][k] = match (&wf_o, &wf_ie) {
                (Some(o), Some(e)) =>  std::cmp::max(o[k-1], e[k-1]) + T::one(),
                (Some(o), None) =>  o[k-1] + T::one(),
                (None, Some(e)) =>  e[k-1] + T::one(),
                _ => wfs_s[1][k],
            }
        }
        for k in bounds[2][0]..bounds[2][1]+1 {
            wfs_s[2][k] = match (&wf_o, &wf_de) {
                (Some(o), Some(e)) => std::cmp::max(o[k+1], e[k+1]),
                (Some(o), None) => o[k+1],
                (None, Some(e)) => e[k+1],
                _ => wfs_s[2][k],
            }
        }
        for k in bounds[0][0]..bounds[0][1]+1 {
            wfs_s[0][k] = match &wf_m {
                Some(m) => *[m[k]+T::one(), wfs_s[1][k], wfs_s[2][k]].iter().max().unwrap(),
                _ => wfs_s[0][k],
            }
        }
    }
    return s;
}

#[allow(dead_code)]
pub unsafe fn wfa_sse4() {
    return;
}

#[allow(dead_code)]
pub unsafe fn wfa_avx2() {
    return;

}