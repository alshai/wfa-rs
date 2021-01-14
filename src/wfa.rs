use num;
use std::option::Option;

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
// T: integer. use larger ints for larger sequences
// q.len() and t.len() must each be smaller than 2^sizeof(T)
// https://doi.org/10.1093/bioinformatics/btaa777
pub fn wfa<T>(opts: WFAOpts, t: &[u8], q: &[u8]) -> usize
    where 
        T: num::PrimInt + num::Signed + num::Zero + num::One + num::FromPrimitive + num::traits::Saturating + std::clone::Clone + std::fmt::Debug,
{
    let max_len = num::pow(2, std::mem::size_of::<T>()*8-1)-1;
    if q.len() > max_len || t.len() > max_len {
        panic!("query length {} or template length {} greater than max: {}. Please consider using a longer integer type", q.len(), t.len(), max_len);
    }
    let k_max = t.len() + q.len() - 2;
    let k_zero = q.len() - 1;
    let a_k = match t.len() >= q.len() { // index of bottom right corner of DP matrix
        true => t.len() - q.len() + k_zero,
        false => q.len() - t.len() + k_zero,
    };
    let a_off : usize = t.len(); // furthest offset of longest diagonal
    // TODO: how can we get around storing <max score> wavefronts?
    let mut mwfs : Vec<Vec<T>> = vec![vec![T::from_i32(-1).unwrap(); k_max+1]; 256]; // match wavefront
    mwfs[0][k_zero] = T::zero();
    let mut iwfs : Vec<Vec<T>> = vec![vec![T::from_i32(-1).unwrap(); k_max+1]; 256]; // insertions
    let mut dwfs : Vec<Vec<T>> = vec![vec![T::from_i32(-1).unwrap(); k_max+1]; 256]; // deletions
    let mut mbounds : [[usize; 2]; 256] = [[k_zero; 2]; 256];
    let mut ibounds : [[usize; 2]; 256] = [[k_zero + 1; 2]; 256];
    let mut dbounds : [[usize; 2]; 256] = [[k_zero - 1; 2]; 256];
    let mut s : usize = 0; // score (must be positive)
    loop {
        // wf_extend(m[s], q, t);
        // extend wavefront by searching down diagonals for matches
        for k in mbounds[s][0]..mbounds[s][1] + 1 {
            if k > k_max {
                break;
            }
            let mut v = (mwfs[s][k].to_i32().unwrap() + (k_zero as i32 - k as i32)) as usize;
            let mut h = mwfs[s][k].to_usize().unwrap_or(0);
            while v < q.len() && h < q.len() && q[v] == t[h] {
                mwfs[s][k] = mwfs[s][k] + T::one();
                v = v + 1;
                h = h + 1;
            }
        }
        // check to see if wavefront extended all the way to bottom right of DP matrix
        if mwfs[s][a_k] >= T::from_usize(a_off).unwrap_or(T::zero()) {
            break;
        }
        // determine if we can start finding gaps
        if s >= (opts.o + opts.e + opts.e) as usize {
            s += *[opts.o + opts.e, opts.x, opts.e].iter().min().unwrap() as usize;
        } else {
            s += *[opts.o + opts.e, opts.x].iter().min().unwrap() as usize;
        }
        // TODO: return error if score goes too long
        // can get rid of this once improvements are made to wavefront storage
        if s > u8::MAX.into() {
            panic!("scores over 255 are currently not supported (score = {})", s);
        }
        // wf_next(m, i, d, q, t, s);
        // get bounds for new wavefront
        // wavefronts for [mismatch, ins, del]
        let mut ixs : [Option<usize>; 3] = [None, None, None];
        // set bounds of M
        if s >= opts.x as usize {
            ixs[0] = Some(s - opts.x as usize);
            mbounds[s][0] = mbounds[ixs[0].unwrap()][0] - 1;
            mbounds[s][1] = mbounds[ixs[0].unwrap()][1] + 1;
        }
        if s >= (opts.o + opts.e) as usize {
            ixs[1] = Some(s - (opts.o + opts.e) as usize);
        }
        // Set bounds of I and D
        if s >= (opts.o + opts.e + opts.e) as usize { // when extension is possible
            ixs[2] = Some(s - opts.e as usize);
            ibounds[s][0] = ibounds[ixs[2].unwrap()][0] - 1;
            dbounds[s][0] = dbounds[ixs[2].unwrap()][0] - 1;
            ibounds[s][1] = ibounds[ixs[2].unwrap()][1] + 1;
            dbounds[s][1] = dbounds[ixs[2].unwrap()][1] + 1;
        }
        // update I, D, M here.
        for k in ibounds[s][0]..ibounds[s][1]+1 {
            iwfs[s][k] = match (ixs[1], ixs[2]) {
                (Some(o), Some(e)) =>  std::cmp::max(mwfs[o][k-1], iwfs[e][k-1]) + T::one(),
                (Some(o), None) =>  mwfs[o][k-1] + T::one(),
                (None, Some(e)) =>  iwfs[e][k-1] + T::one(),
                _ => iwfs[s][k],
            }
        }
        for k in dbounds[s][0]..dbounds[s][1]+1 {
            dwfs[s][k] = match (ixs[1], ixs[2]) {
                (Some(o), Some(e)) => std::cmp::max(mwfs[o][k+1], dwfs[e][k+1]),
                (Some(o), None) => mwfs[o][k+1],
                (None, Some(e)) => dwfs[e][k+1],
                _ => dwfs[s][k],
            }
        }
        for k in mbounds[s][0]..mbounds[s][1]+1 {
            mwfs[s][k] = match ixs[0] {
                Some(i) => *[mwfs[i][k]+T::one(), iwfs[s][k], dwfs[s][k]].iter().max().unwrap(),
                _ => mwfs[s][k],
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