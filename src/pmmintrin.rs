// Copyright 2016 Eli Friedman.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#![allow(unused_variables)]

use conversions::Convert128;
use simd::i8x16;
use __m128;
use __m128i;
use __m128d;

// Declarations copied from the llvmint crate.
#[allow(improper_ctypes)]
extern {
    #[link_name = "llvm.x86.sse3.ldu.dq"]
    pub fn sse3_ldu_dq(a: *mut i8) -> i8x16;
}

/// addsubpd
#[inline]
pub fn _mm_addsub_pd (a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// addsubps
#[inline]
pub fn _mm_addsub_ps (a: __m128, b: __m128) -> __m128 {
    unimplemented!()
}
/// haddpd
#[inline]
pub fn _mm_hadd_pd (a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// haddps
#[inline]
pub fn _mm_hadd_ps (a: __m128, b: __m128) -> __m128 {
    unimplemented!()
}
/// hsubpd
#[inline]
pub fn _mm_hsub_pd (a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// hsubps
#[inline]
pub fn _mm_hsub_ps (a: __m128, b: __m128) -> __m128 {
    unimplemented!()
}
/// lddqu
#[inline]
pub unsafe fn _mm_lddqu_si128(mem_addr: *const __m128i) -> __m128i {
    sse3_ldu_dq(mem_addr as *mut i8).as_i64x2()
}
/// movddup
#[inline]
pub unsafe fn _mm_loaddup_pd (mem_addr: *const f64) ->  __m128d {
    unimplemented!()
}
/// movddup
#[inline]
pub fn _mm_movedup_pd (a: __m128d) ->  __m128d {
    __m128d::new(a.extract(0), a.extract(0))
}
/// movshdup
#[inline]
pub fn _mm_movehdup_ps (a: __m128) -> __m128 {
    unimplemented!()
}
/// movsldup
#[inline]
pub fn _mm_moveldup_ps (a: __m128) -> __m128 {
    unimplemented!()
}
