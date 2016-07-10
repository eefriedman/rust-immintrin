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
use __m128d;
use simd::i32x4;
use simd::i16x8;
use simd::i8x16;
use simd::u8x16;

pub trait Convert128 : Sized {
    fn as_f32x4(self) -> __m128;
    fn as_i64x2(self) -> __m128i;
    fn as_f64x2(self) -> __m128d;
    fn as_i32x4(self) -> i32x4;
    fn as_i16x8(self) -> i16x8;
    fn as_i8x16(self) -> i8x16;
    fn as_u8x16(self) -> u8x16;
}

macro_rules! impl_convert128 {
    ( $x:ident ) => {
        impl Convert128 for $x {
            fn as_f32x4(self) -> __m128 { unsafe { transmute(self) } }
            fn as_i64x2(self) -> __m128i { unsafe { transmute(self) } }
            fn as_f64x2(self) -> __m128d { unsafe { transmute(self) } }
            fn as_i32x4(self) -> i32x4 { unsafe { transmute(self) } }
            fn as_i16x8(self) -> i16x8 { unsafe { transmute(self) } }
            fn as_i8x16(self) -> i8x16 { unsafe { transmute(self) } }
            fn as_u8x16(self) -> u8x16 { unsafe { transmute(self) } }
        }
    };
}

impl_convert128!(__m128);
impl_convert128!(__m128i);
impl_convert128!(__m128d);
impl_convert128!(i32x4);
impl_convert128!(i16x8);
impl_convert128!(i8x16);
