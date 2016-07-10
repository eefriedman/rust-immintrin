// Copyright 2016 Eli Friedman.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#![allow(unused_variables)]

use conversions::Convert128;
use simd::x86::ssse3::{Ssse3U8x16, Ssse3I8x16, Ssse3I16x8, Ssse3I32x4};
use __m128i;

/// pabsw
#[inline]
pub fn _mm_abs_epi16(a: __m128i) -> __m128i {
    a.as_i16x8().abs().as_i64x2()
}
/// pabsd
#[inline]
pub fn _mm_abs_epi32(a: __m128i) -> __m128i {
    a.as_i32x4().abs().as_i64x2()
}
/// pabsb
#[inline]
pub fn _mm_abs_epi8(a: __m128i) -> __m128i {
    a.as_i8x16().abs().as_i64x2()
}
/// palignr
#[inline]
pub fn _mm_alignr_epi8(a: __m128i, b: __m128i, count: i32) -> __m128i {
    unimplemented!()
}
/// phaddw
#[inline]
pub fn _mm_hadd_epi16(a: __m128i, b: __m128i) -> __m128i {
    a.as_i16x8().hadd(b.as_i16x8()).as_i64x2()
}
/// phaddd
#[inline]
pub fn _mm_hadd_epi32(a: __m128i, b: __m128i) -> __m128i {
    a.as_i32x4().hadd(b.as_i32x4()).as_i64x2()
}
/// phaddsw
#[inline]
pub fn _mm_hadds_epi16(a: __m128i, b: __m128i) -> __m128i {
    a.as_i16x8().hadds(b.as_i16x8()).as_i64x2()
}
/// phsubw
#[inline]
pub fn _mm_hsub_epi16(a: __m128i, b: __m128i) -> __m128i {
    a.as_i16x8().hsub(b.as_i16x8()).as_i64x2()
}
/// phsubd
#[inline]
pub fn _mm_hsub_epi32(a: __m128i, b: __m128i) -> __m128i {
    a.as_i32x4().hsub(b.as_i32x4()).as_i64x2()
}
/// phsubsw
#[inline]
pub fn _mm_hsubs_epi16(a: __m128i, b: __m128i) -> __m128i {
    a.as_i16x8().hsubs(b.as_i16x8()).as_i64x2()
}
/// pmaddubsw
#[inline]
pub fn _mm_maddubs_epi16(a: __m128i, b: __m128i) -> __m128i {
    // FIXME: Missing inline marking in SIMD crate
    a.as_u8x16().maddubs(b.as_i8x16()).as_i64x2()
}
/// pmulhrsw
#[inline]
pub fn _mm_mulhrs_epi16(a: __m128i, b: __m128i) -> __m128i {
    a.as_i16x8().mulhrs(b.as_i16x8()).as_i64x2()
}
/// pshufb
#[inline]
pub fn _mm_shuffle_epi8(a: __m128i, b: __m128i) -> __m128i {
    a.as_i8x16().shuffle_bytes(b.as_i8x16()).as_i64x2()
}
/// psignw
#[inline]
pub fn _mm_sign_epi16(a: __m128i, b: __m128i) -> __m128i {
    a.as_i16x8().sign(b.as_i16x8()).as_i64x2()
}
/// psignd
#[inline]
pub fn _mm_sign_epi32(a: __m128i, b: __m128i) -> __m128i {
    a.as_i32x4().sign(b.as_i32x4()).as_i64x2()
}
/// psignb
#[inline]
pub fn _mm_sign_epi8(a: __m128i, b: __m128i) -> __m128i {
    a.as_i8x16().sign(b.as_i8x16()).as_i64x2()
}

/// The methods in this module can't be implemented because Rust doesn't
/// expose the LLVM x86_mmx type.
pub mod unimplemented_mmx {
    use __m64;
    /// pabsw
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_abs_pi16(a: __m64) -> __m64 {
        unimplemented!()
    }
    /// pabsd
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_abs_pi32(a: __m64) -> __m64 {
        unimplemented!()
    }
    /// pabsb
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_abs_pi8(a: __m64) -> __m64 {
        unimplemented!()
    }
    /// palignr
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_alignr_pi8(a: __m64, b: __m64, count: i32) -> __m64 {
        unimplemented!()
    }
    /// phaddw
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_hadd_pi16(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// phaddw
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_hadd_pi32(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// phaddsw
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_hadds_pi16(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// phsubw
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_hsub_pi16(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// phsubd
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_hsub_pi32(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// phsubswb: __m64
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_hsubs_pi16(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// pmaddubsw
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_maddubs_pi16(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// pmulhrsw
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_mulhrs_pi16(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// pshufb
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_shuffle_pi8(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// psignw
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_sign_pi16(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// psignd
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_sign_pi32(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// psignb
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_sign_pi8(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
}