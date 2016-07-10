// Copyright 2016 Eli Friedman.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use std::mem::transmute;
use __m128;
use __m128i;

#[inline]
pub fn _mm_castsi128_ps(a: __m128i) -> __m128 {
    unsafe { transmute(a) }
}

#[inline]
pub fn _mm_castps_si128(a: __m128) -> __m128i {
    unsafe { transmute(a) }
}