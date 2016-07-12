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
use simd::Simd;

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

// Copied from SIMD crate.
#[allow(dead_code)]
extern "platform-intrinsic" {
    fn simd_eq<T: Simd<Bool = U>, U>(x: T, y: T) -> U;
    fn simd_ne<T: Simd<Bool = U>, U>(x: T, y: T) -> U;
    fn simd_lt<T: Simd<Bool = U>, U>(x: T, y: T) -> U;
    fn simd_le<T: Simd<Bool = U>, U>(x: T, y: T) -> U;
    fn simd_gt<T: Simd<Bool = U>, U>(x: T, y: T) -> U;
    fn simd_ge<T: Simd<Bool = U>, U>(x: T, y: T) -> U;

    fn simd_shuffle2<T: Simd, U: Simd<Elem = T::Elem>>(x: T, y: T, idx: [u32; 2]) -> U;
    fn simd_shuffle4<T: Simd, U: Simd<Elem = T::Elem>>(x: T, y: T, idx: [u32; 4]) -> U;
    fn simd_shuffle8<T: Simd, U: Simd<Elem = T::Elem>>(x: T, y: T, idx: [u32; 8]) -> U;
    fn simd_shuffle16<T: Simd, U: Simd<Elem = T::Elem>>(x: T, y: T, idx: [u32; 16]) -> U;

    fn simd_insert<T: Simd<Elem = U>, U>(x: T, idx: u32, val: U) -> T;
    fn simd_extract<T: Simd<Elem = U>, U>(x: T, idx: u32) -> U;

    fn simd_cast<T: Simd, U: Simd>(x: T) -> U;

    fn simd_add<T: Simd>(x: T, y: T) -> T;
    fn simd_sub<T: Simd>(x: T, y: T) -> T;
    fn simd_mul<T: Simd>(x: T, y: T) -> T;
    fn simd_div<T: Simd>(x: T, y: T) -> T;
    fn simd_shl<T: Simd>(x: T, y: T) -> T;
    fn simd_shr<T: Simd>(x: T, y: T) -> T;
    fn simd_and<T: Simd>(x: T, y: T) -> T;
    fn simd_or<T: Simd>(x: T, y: T) -> T;
    fn simd_xor<T: Simd>(x: T, y: T) -> T;
}