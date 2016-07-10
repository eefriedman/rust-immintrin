// Copyright 2016 Eli Friedman.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#![feature(repr_simd)]
#![feature(link_llvm_intrinsics)]
#![feature(platform_intrinsics)]
#![feature(simd_ffi)]
#![feature(cfg_target_feature)]

#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

extern crate simd;

#[cfg(target_feature="sse2")]
use simd::f32x4;
#[cfg(target_feature="sse2")]
use simd::x86::sse2::{f64x2, i64x2};

#[cfg(target_feature="sse2")]
pub type __m128 = f32x4;
#[cfg(target_feature="sse2")]
pub type __m128d = f64x2;
#[cfg(target_feature="sse2")]
pub type __m128i = i64x2;
#[repr(simd)]
#[derive(Debug, Copy, Clone)]
pub struct __m64(i64);

#[cfg(target_feature="sse2")]
pub mod xmmintrin;
#[cfg(target_feature="sse2")]
pub mod emmintrin;
#[cfg(target_feature="sse3")]
pub mod pmmintrin;
#[cfg(target_feature="ssse3")]
pub mod tmmintrin;
#[cfg(target_feature="sse2")]
mod conversions;
