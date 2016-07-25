// Copyright 2016 Eli Friedman.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#![allow(unused_variables)]

use conversions::Convert128;
use simd::x86::sse2::{Sse2I8x16, Sse2U8x16, Sse2I16x8, Sse2U16x8, Sse2F32x4, Sse2F64x2};
use __m128;
use __m128i;
use __m128d;
use simd::i32x4;
use simd::i16x8;
use simd::x86::sse2::bool64fx2;
use std::ptr::copy_nonoverlapping;
use std::mem::transmute;
use simd_shuffle4;
use simd_shuffle8;

// Declarations copied from the llvmint crate.
#[allow(improper_ctypes)]
extern {
    #[link_name = "llvm.x86.sse2.sqrt.sd"]
    pub fn sse2_sqrt_sd(a: __m128d) -> __m128d;
    #[link_name = "llvm.x86.sse2.cvtsd2si"]
    pub fn sse2_cvtsd2si(a: __m128d) -> i32;
    #[link_name = "llvm.x86.sse2.cvttsd2si"]
    pub fn sse2_cvttsd2si(a: __m128d) -> i32;
    #[link_name = "llvm.x86.sse2.cvtsd2si64"]
    pub fn sse2_cvtsd2si64(a: __m128d) -> i64;
    #[link_name = "llvm.x86.sse2.cvttsd2si64"]
    pub fn sse2_cvttsd2si64(a: __m128d) -> i64;
    #[link_name = "llvm.x86.sse2.cvtpd2ps"]
    pub fn sse2_cvtpd2ps(a: __m128d) -> __m128;
    #[link_name = "llvm.x86.sse2.cvtps2dq"]
    pub fn sse2_cvtps2dq(a: __m128) -> i32x4;
    #[link_name = "llvm.x86.sse2.cvttps2dq"]
    pub fn sse2_cvttps2dq(a: __m128) -> i32x4;
    #[link_name = "llvm.x86.sse2.mfence"]
    pub fn sse2_mfence() ->();
    #[link_name = "llvm.x86.sse2.min.sd"]
    pub fn sse2_min_sd(a: __m128d, b: __m128d) -> __m128d;
    #[link_name = "llvm.x86.sse2.max.sd"]
    pub fn sse2_max_sd(a: __m128d, b: __m128d) -> __m128d;
    #[link_name = "llvm.x86.sse2.cmp.sd"]
    pub fn sse2_cmp_sd(a: __m128d, b: __m128d, c: i8) -> __m128d;
    #[link_name = "llvm.x86.sse2.comieq.sd"]
    pub fn sse2_comieq_sd(a: __m128d, b: __m128d) -> i32;
    #[link_name = "llvm.x86.sse2.comilt.sd"]
    pub fn sse2_comilt_sd(a: __m128d, b: __m128d) -> i32;
    #[link_name = "llvm.x86.sse2.comile.sd"]
    pub fn sse2_comile_sd(a: __m128d, b: __m128d) -> i32;
    #[link_name = "llvm.x86.sse2.comigt.sd"]
    pub fn sse2_comigt_sd(a: __m128d, b: __m128d) -> i32;
    #[link_name = "llvm.x86.sse2.comige.sd"]
    pub fn sse2_comige_sd(a: __m128d, b: __m128d) -> i32;
    #[link_name = "llvm.x86.sse2.comineq.sd"]
    pub fn sse2_comineq_sd(a: __m128d, b: __m128d) -> i32;
    #[link_name = "llvm.x86.sse2.ucomieq.sd"]
    pub fn sse2_ucomieq_sd(a: __m128d, b: __m128d) -> i32;
    #[link_name = "llvm.x86.sse2.ucomilt.sd"]
    pub fn sse2_ucomilt_sd(a: __m128d, b: __m128d) -> i32;
    #[link_name = "llvm.x86.sse2.ucomile.sd"]
    pub fn sse2_ucomile_sd(a: __m128d, b: __m128d) -> i32;
    #[link_name = "llvm.x86.sse2.ucomigt.sd"]
    pub fn sse2_ucomigt_sd(a: __m128d, b: __m128d) -> i32;
    #[link_name = "llvm.x86.sse2.ucomige.sd"]
    pub fn sse2_ucomige_sd(a: __m128d, b: __m128d) -> i32;
    #[link_name = "llvm.x86.sse2.ucomineq.sd"]
    pub fn sse2_ucomineq_sd(a: __m128d, b: __m128d) -> i32;
    #[link_name = "llvm.x86.sse2.cvtsi2sd"]
    pub fn sse2_cvtsi2sd(a: __m128d, b: i32) -> __m128d;
    #[link_name = "llvm.x86.sse2.cvtsi642sd"]
    pub fn sse2_cvtsi642sd(a: __m128d, b: i64) -> __m128d;
    #[link_name = "llvm.x86.sse2.pslli.w"]
    pub fn sse2_pslli_w(a: i16x8, b: i32) -> i16x8;
    #[link_name = "llvm.x86.sse2.pslli.d"]
    pub fn sse2_pslli_d(a: i32x4, b: i32) -> i32x4;
    #[link_name = "llvm.x86.sse2.pslli.q"]
    pub fn sse2_pslli_q(a: __m128i, b: i32) -> __m128i;
    #[link_name = "llvm.x86.sse2.psrli.w"]
    pub fn sse2_psrli_w(a: i16x8, b: i32) -> i16x8;
    #[link_name = "llvm.x86.sse2.psrli.d"]
    pub fn sse2_psrli_d(a: i32x4, b: i32) -> i32x4;
    #[link_name = "llvm.x86.sse2.psrli.q"]
    pub fn sse2_psrli_q(a: __m128i, b: i32) -> __m128i;
    #[link_name = "llvm.x86.sse2.psrai.w"]
    pub fn sse2_psrai_w(a: i16x8, b: i32) -> i16x8;
    #[link_name = "llvm.x86.sse2.psrai.d"]
    pub fn sse2_psrai_d(a: i32x4, b: i32) -> i32x4;
    #[link_name = "llvm.x86.sse2.psll.dq"]
    pub fn sse2_psll_dq(a: __m128i, b: i32) -> __m128i;
    #[link_name = "llvm.x86.sse2.psrl.dq"]
    pub fn sse2_psrl_dq(a: __m128i, b: i32) -> __m128i;
    #[link_name = "llvm.x86.sse2.psll.w"]
    pub fn sse2_psll_w(a: i16x8, b: i16x8) -> i16x8;
    #[link_name = "llvm.x86.sse2.psll.d"]
    pub fn sse2_psll_d(a: i32x4, b: i32x4) -> i32x4;
    #[link_name = "llvm.x86.sse2.psll.q"]
    pub fn sse2_psll_q(a: __m128i, b: __m128i) -> __m128i;
    #[link_name = "llvm.x86.sse2.psrl.w"]
    pub fn sse2_psrl_w(a: i16x8, b: i16x8) -> i16x8;
    #[link_name = "llvm.x86.sse2.psrl.d"]
    pub fn sse2_psrl_d(a: i32x4, b: i32x4) -> i32x4;
    #[link_name = "llvm.x86.sse2.psrl.q"]
    pub fn sse2_psrl_q(a: __m128i, b: __m128i) -> __m128i;
    #[link_name = "llvm.x86.sse2.psra.w"]
    pub fn sse2_psra_w(a: i16x8, b: i16x8) -> i16x8;
    #[link_name = "llvm.x86.sse2.psra.d"]
    pub fn sse2_psra_d(a: i32x4, b: i32x4) -> i32x4;
    #[link_name = "llvm.x86.sse2.clflush"]
    pub fn sse2_clflush(a: *mut i8) -> ();
    #[link_name = "llvm.x86.sse2.cvtdq2pd"]
    pub fn sse2_cvtdq2pd(a: i32x4) -> __m128d;
    #[link_name = "llvm.x86.sse2.cvtpd2dq"]
    pub fn sse2_cvtpd2dq(a: __m128d) -> i32x4;
    #[link_name = "llvm.x86.sse2.cvttpd2dq"]
    pub fn sse2_cvttpd2dq(a: __m128d) -> i32x4;
    #[link_name = "llvm.x86.sse2.cvtsd2ss"]
    pub fn sse2_cvtsd2ss(a: __m128, b: __m128d) -> __m128;
    #[link_name = "llvm.x86.sse2.cvtss2sd"]
    pub fn sse2_cvtss2sd(a: __m128d, b: __m128) -> __m128d;
}

fn convert_bool64fx2_to_m128d(a: bool64fx2) -> __m128d {
    unsafe { transmute(a) }
}

/// paddw
#[inline]
pub fn _mm_add_epi16(a: __m128i, b: __m128i) -> __m128i {
    (a.as_i16x8() + b.as_i16x8()).as_i64x2()
}
/// paddd
#[inline]
pub fn _mm_add_epi32(a: __m128i, b: __m128i) -> __m128i {
    (a.as_i32x4() + b.as_i32x4()).as_i64x2()
}
/// paddq
#[inline]
pub fn _mm_add_epi64(a: __m128i, b: __m128i) -> __m128i {
    a + b
}
/// paddb
#[inline]
pub fn _mm_add_epi8(a: __m128i, b: __m128i) -> __m128i {
    (a.as_i8x16() + b.as_i8x16()).as_i64x2()
}
/// addpd
#[inline]
pub fn _mm_add_pd(a: __m128d, b: __m128d) -> __m128d {
    a + b
}
/// addsd
#[inline]
pub fn _mm_add_sd(a: __m128d, b: __m128d) -> __m128d {
    a.replace(0, a.extract(0) + b.extract(0))
}
/// paddsw
#[inline]
pub fn _mm_adds_epi16(a: __m128i, b: __m128i) -> __m128i {
    a.as_i16x8().adds(b.as_i16x8()).as_i64x2()
}
/// paddsb
#[inline]
pub fn _mm_adds_epi8(a: __m128i, b: __m128i) -> __m128i {
    a.as_i8x16().adds(b.as_i8x16()).as_i64x2()
}
/// paddusw
#[inline]
pub fn _mm_adds_epu16(a: __m128i, b: __m128i) -> __m128i {
    a.as_u16x8().adds(b.as_u16x8()).as_i64x2()
}
/// paddusb
#[inline]
pub fn _mm_adds_epu8(a: __m128i, b: __m128i) -> __m128i {
    a.as_u8x16().adds(b.as_u8x16()).as_i64x2()
}
/// andpd
#[inline]
pub fn _mm_and_pd(a: __m128d, b: __m128d) -> __m128d {
    (a.as_i64x2() & b.as_i64x2()).as_f64x2()
}
/// pand
#[inline]
pub fn _mm_and_si128(a: __m128i, b: __m128i) -> __m128i {
    a & b
}
/// andnpd
#[inline]
pub fn _mm_andnot_pd(a: __m128d, b: __m128d) -> __m128d {
    // FIXME: Not operator from simd doesn't get inlined?!
    ((a.as_i64x2() ^ __m128i::splat(!0)) & b.as_i64x2()).as_f64x2()
}
/// pandn
#[inline]
pub fn _mm_andnot_si128(a: __m128i, b: __m128i) -> __m128i {
    // FIXME: Not operator from simd doesn't get inlined?!
    (a ^ __m128i::splat(!0)) & b
}
/// pavgw
#[inline]
pub fn _mm_avg_epu16(a: __m128i, b: __m128i) -> __m128i {
    a.as_u16x8().avg(b.as_u16x8()).as_i64x2()
}
/// pavgb
#[inline]
pub fn _mm_avg_epu8(a: __m128i, b: __m128i) -> __m128i {
    a.as_u8x16().avg(b.as_u8x16()).as_i64x2()
}
/// pslldq
#[inline]
pub fn _mm_bslli_si128(a: __m128i, imm8: i32) -> __m128i {
    _mm_slli_si128(a, imm8)
}
/// psrldq
#[inline]
pub fn _mm_bsrli_si128(a: __m128i, imm8: i32) -> __m128i {
    _mm_srli_si128(a, imm8)
}
#[inline]
pub fn _mm_castpd_ps(a: __m128d) -> __m128 {
    a.as_f32x4()
}
#[inline]
pub fn _mm_castpd_si128(a: __m128d) -> __m128i {
    a.as_i64x2()
}
#[inline]
pub fn _mm_castps_pd(a: __m128) -> __m128d {
    a.as_f64x2()
}
#[inline]
pub fn _mm_castps_si128(a: __m128) -> __m128i {
    a.as_i64x2()
}
#[inline]
pub fn _mm_castsi128_pd(a: __m128i) -> __m128d {
    a.as_f64x2()
}
#[inline]
pub fn _mm_castsi128_ps(a: __m128i) -> __m128 {
    a.as_f32x4()
}
/// clflush
#[inline]
pub fn _mm_clflush(p: *const u8) {
    unsafe { sse2_clflush(p as *mut i8) }
}
/// pcmpeqw
#[inline]
pub fn _mm_cmpeq_epi16(a: __m128i, b: __m128i) -> __m128i {
    a.as_i16x8().eq(b.as_i16x8()).to_repr().as_i64x2()
}
/// pcmpeqd
#[inline]
pub fn _mm_cmpeq_epi32(a: __m128i, b: __m128i) -> __m128i {
    a.as_i32x4().eq(b.as_i32x4()).to_repr().as_i64x2()
}
/// pcmpeqb
#[inline]
pub fn _mm_cmpeq_epi8(a: __m128i, b: __m128i) -> __m128i {
    a.as_i8x16().eq(b.as_i8x16()).to_repr().as_i64x2()
}
/// cmppd
#[inline]
pub fn _mm_cmpeq_pd(a: __m128d, b: __m128d) -> __m128d {
    convert_bool64fx2_to_m128d(a.eq(b))
}
/// cmpsd
#[inline]
pub fn _mm_cmpeq_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpge_pd(a: __m128d, b: __m128d) -> __m128d {
    convert_bool64fx2_to_m128d(a.ge(b))
}
/// cmpsd
#[inline]
pub fn _mm_cmpge_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// pcmpgtw
#[inline]
pub fn _mm_cmpgt_epi16(a: __m128i, b: __m128i) -> __m128i {
    a.as_i16x8().gt(b.as_i16x8()).to_repr().as_i64x2()
}
/// pcmpgtd
#[inline]
pub fn _mm_cmpgt_epi32(a: __m128i, b: __m128i) -> __m128i {
    a.as_i32x4().gt(b.as_i32x4()).to_repr().as_i64x2()
}
/// pcmpgtb
#[inline]
pub fn _mm_cmpgt_epi8(a: __m128i, b: __m128i) -> __m128i {
    a.as_i8x16().gt(b.as_i8x16()).to_repr().as_i64x2()
}
/// cmppd
#[inline]
pub fn _mm_cmpgt_pd(a: __m128d, b: __m128d) -> __m128d {
    convert_bool64fx2_to_m128d(a.gt(b))
}
/// cmpsd
#[inline]
pub fn _mm_cmpgt_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmple_pd(a: __m128d, b: __m128d) -> __m128d {
    convert_bool64fx2_to_m128d(a.le(b))
}
/// cmpsd
#[inline]
pub fn _mm_cmple_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// pcmpgtw
#[inline]
pub fn _mm_cmplt_epi16(a: __m128i, b: __m128i) -> __m128i {
    _mm_cmpgt_epi16(b, a)
}
/// pcmpgtd
#[inline]
pub fn _mm_cmplt_epi32(a: __m128i, b: __m128i) -> __m128i {
    _mm_cmpgt_epi32(b, a)
}
/// pcmpgtb
#[inline]
pub fn _mm_cmplt_epi8(a: __m128i, b: __m128i) -> __m128i {
    _mm_cmpgt_epi8(b, a)
}
/// cmppd
#[inline]
pub fn _mm_cmplt_pd(a: __m128d, b: __m128d) -> __m128d {
    convert_bool64fx2_to_m128d(a.lt(b))
}
/// cmpsd
#[inline]
pub fn _mm_cmplt_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpneq_pd(a: __m128d, b: __m128d) -> __m128d {
    convert_bool64fx2_to_m128d(a.ne(b))
}
/// cmpsd
#[inline]
pub fn _mm_cmpneq_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpnge_pd(a: __m128d, b: __m128d) -> __m128d {
    convert_bool64fx2_to_m128d(a.lt(b) | (a.ne(a) | b.ne(b)))
}
/// cmpsd
#[inline]
pub fn _mm_cmpnge_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpngt_pd(a: __m128d, b: __m128d) -> __m128d {
    convert_bool64fx2_to_m128d(a.le(b) | (a.ne(a) | b.ne(b)))
}
/// cmpsd
#[inline]
pub fn _mm_cmpngt_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpnle_pd(a: __m128d, b: __m128d) -> __m128d {
    convert_bool64fx2_to_m128d(a.gt(b) | (a.ne(a) | b.ne(b)))
}
/// cmpsd
#[inline]
pub fn _mm_cmpnle_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpnlt_pd(a: __m128d, b: __m128d) -> __m128d {
    convert_bool64fx2_to_m128d(a.ge(b) | (a.ne(a) | b.ne(b)))
}
/// cmpsd
#[inline]
pub fn _mm_cmpnlt_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpord_pd(a: __m128d, b: __m128d) -> __m128d {
    convert_bool64fx2_to_m128d(a.eq(a) & b.eq(b))
}
/// cmpsd
#[inline]
pub fn _mm_cmpord_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpunord_pd(a: __m128d, b: __m128d) -> __m128d {
    convert_bool64fx2_to_m128d(a.ne(a) | b.ne(b))
}
/// cmpsd
#[inline]
pub fn _mm_cmpunord_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// comisd
#[inline]
pub fn _mm_comieq_sd(a: __m128d, b: __m128d) -> i32 {
    unsafe { sse2_comieq_sd(a, b) }
}
/// comisd
#[inline]
pub fn _mm_comige_sd(a: __m128d, b: __m128d) -> i32 {
    unsafe { sse2_comige_sd(a, b) }
}
/// comisd
#[inline]
pub fn _mm_comigt_sd(a: __m128d, b: __m128d) -> i32 {
    unsafe { sse2_comigt_sd(a, b) }
}
/// comisd
#[inline]
pub fn _mm_comile_sd(a: __m128d, b: __m128d) -> i32 {
    unsafe { sse2_comile_sd(a, b) }
}
/// comisd
#[inline]
pub fn _mm_comilt_sd(a: __m128d, b: __m128d) -> i32 {
    unsafe { sse2_comilt_sd(a, b) }
}
/// comisd
#[inline]
pub fn _mm_comineq_sd(a: __m128d, b: __m128d) -> i32 {
    unsafe { sse2_comineq_sd(a, b) }
}
/// cvtdq2pd
#[inline]
pub fn _mm_cvtepi32_pd(a: __m128i) -> __m128d {
    unsafe { sse2_cvtdq2pd(a.as_i32x4()) }
}
/// cvtdq2ps
#[inline]
pub fn _mm_cvtepi32_ps(a: __m128i) -> __m128 {
    a.as_i32x4().to_f32()
}
/// cvtpd2dq
#[inline]
pub fn _mm_cvtpd_epi32(a: __m128d) -> __m128i {
    unsafe { sse2_cvtpd2dq(a).as_i64x2() }
}
/// cvtpd2ps
#[inline]
pub fn _mm_cvtpd_ps(a: __m128d) -> __m128 {
    unsafe { sse2_cvtpd2ps(a) }
}
/// cvtps2dq
#[inline]
pub fn _mm_cvtps_epi32(a: __m128) -> __m128i {
    unsafe { sse2_cvtps2dq(a).as_i64x2() }
}
/// cvtps2pd
#[inline]
pub fn _mm_cvtps_pd(a: __m128) -> __m128d {
    a.to_f64()
}
/// movsd
#[inline]
pub fn _mm_cvtsd_f64(a: __m128d) -> f64 {
    a.extract(0)
}
/// cvtsd2si
#[inline]
pub fn _mm_cvtsd_si32(a: __m128d) -> i32 {
    unsafe { sse2_cvtsd2si(a) }
}
/// cvtsd2si
#[inline]
#[cfg(target_pointer_width = "64")]
pub fn _mm_cvtsd_si64(a: __m128d) -> i64 {
    unsafe { sse2_cvtsd2si64(a) }
}
/// cvtsd2si
#[inline]
#[cfg(target_pointer_width = "64")]
pub fn _mm_cvtsd_si64x(a: __m128d) -> i64 {
    _mm_cvtsd_si64(a)
}
/// cvtsd2ss
#[inline]
pub fn _mm_cvtsd_ss(a: __m128, b: __m128d) -> __m128 {
    unsafe { sse2_cvtsd2ss(a, b) }
}
/// movd
#[inline]
pub fn _mm_cvtsi128_si32(a: __m128i) -> i32 {
    a.as_i32x4().extract(0)
}
/// movq
#[inline]
pub fn _mm_cvtsi128_si64(a: __m128i) -> i64 {
    a.extract(0)
}
/// movq
#[inline]
pub fn _mm_cvtsi128_si64x(a: __m128i) -> i64 {
    _mm_cvtsi128_si64(a)
}
/// cvtsi2sd
#[inline]
pub fn _mm_cvtsi32_sd(a: __m128d, b: i32) -> __m128d {
    unsafe { sse2_cvtsi2sd(a, b) }
}
/// movd
#[inline]
pub fn _mm_cvtsi32_si128(a: i32) -> __m128i {
    i32x4::new(a, 0, 0, 0).as_i64x2()
}
/// cvtsi2sd
#[inline]
#[cfg(target_pointer_width = "64")]
pub fn _mm_cvtsi64_sd(a: __m128d, b: i64) -> __m128d {
    unsafe { sse2_cvtsi642sd(a, b) }
}
/// movq
#[inline]
pub fn _mm_cvtsi64_si128(a: i64) -> __m128i {
    __m128i::new(a, 0)
}
/// cvtsi2sd
#[inline]
#[cfg(target_pointer_width = "64")]
pub fn _mm_cvtsi64x_sd(a: __m128d, b: i64) -> __m128d {
    _mm_cvtsi64_sd(a, b)
}
/// movq
#[inline]
pub fn _mm_cvtsi64x_si128(a: i64) -> __m128i {
    _mm_cvtsi64_si128(a)
}
/// cvtss2sd
#[inline]
pub fn _mm_cvtss_sd(a: __m128d, b: __m128) -> __m128d {
    unsafe { sse2_cvtss2sd(a, b) }
}
/// cvttpd2dq
#[inline]
pub fn _mm_cvttpd_epi32(a: __m128d) -> __m128i {
    unsafe { sse2_cvttpd2dq(a).as_i64x2() }
}
/// cvttps2dq
#[inline]
pub fn _mm_cvttps_epi32(a: __m128) -> __m128i {
    unsafe { sse2_cvttps2dq(a).as_i64x2() }
}
/// cvttsd2si
#[inline]
pub fn _mm_cvttsd_si32(a: __m128d) -> i32 {
    unsafe { sse2_cvttsd2si(a) }
}
/// cvttsd2si
#[inline]
#[cfg(target_pointer_width = "64")]
pub fn _mm_cvttsd_si64(a: __m128d) -> i64 {
    unsafe { sse2_cvttsd2si64(a) }
}
/// cvttsd2si
#[inline]
#[cfg(target_pointer_width = "64")]
pub fn _mm_cvttsd_si64x(a: __m128d) -> i64 {
    _mm_cvttsd_si64(a)
}
/// divpd
#[inline]
pub fn _mm_div_pd(a: __m128d, b: __m128d) -> __m128d {
    a / b
}
/// divsd
#[inline]
pub fn _mm_div_sd(a: __m128d, b: __m128d) -> __m128d {
    a.replace(0, a.extract(0) / b.extract(0))
}
/// pextrw
#[inline]
pub fn _mm_extract_epi16(a: __m128i, imm8: i32) -> i32 {
    a.as_i16x8().extract(imm8 as u32) as i32
}
/// pinsrw
#[inline]
pub fn _mm_insert_epi16(a: __m128i, i: i32, imm8: i32) -> __m128i {
    a.as_i16x8().replace(imm8 as u32, i as i16).as_i64x2()
}
/// lfence
#[inline]
pub fn _mm_lfence() {
    unimplemented!()
}
/// movapd
#[inline]
pub unsafe fn _mm_load_pd(mem_addr: *const f64) -> __m128d {
    unimplemented!()
}
#[inline]
pub unsafe fn _mm_load_pd1(mem_addr: *const f64) -> __m128d {
    unimplemented!()
}
/// movsd
#[inline]
pub unsafe fn _mm_load_sd(mem_addr: *const f64) -> __m128d {
    unimplemented!()
}
/// movdqa
#[inline]
pub unsafe fn _mm_load_si128(mem_addr: *const __m128i) -> __m128i {
    unimplemented!()
}
#[inline]
pub unsafe fn _mm_load1_pd(mem_addr: *const f64) -> __m128d {
    unimplemented!()
}
/// movhpd
#[inline]
pub unsafe fn _mm_loadh_pd(a: __m128d, mem_addr: *const f64) -> __m128d {
    unimplemented!()
}
/// movq
#[inline]
pub unsafe fn _mm_loadl_epi64(mem_addr: *const __m128i) -> __m128i {
    unimplemented!()
}
/// movlpd
#[inline]
pub unsafe fn _mm_loadl_pd(a: __m128d, mem_addr: *const f64) -> __m128d {
    unimplemented!()
}
#[inline]
pub unsafe fn _mm_loadr_pd(mem_addr: *const f64) -> __m128d {
    unimplemented!()
}
/// movupd
#[inline]
pub unsafe fn _mm_loadu_pd(mem_addr: *const f64) -> __m128d {
    unimplemented!()
}
/// movdqu
#[inline]
pub unsafe fn _mm_loadu_si128(mem_addr: *const __m128i) -> __m128i {
    unimplemented!()
}
/// pmaddwd
#[inline]
pub fn _mm_madd_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// maskmovdqu
#[inline]
pub unsafe fn _mm_maskmoveu_si128(a: __m128i, mask: __m128i, mem_addr: *mut i8) {
    unimplemented!()
}
/// pmaxsw
#[inline]
pub fn _mm_max_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// pmaxub
#[inline]
pub fn _mm_max_epu8(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// maxpd
#[inline]
pub fn _mm_max_pd(a: __m128d, b: __m128d) -> __m128d {
    a.max(b)
}
/// maxsd
#[inline]
pub fn _mm_max_sd(a: __m128d, b: __m128d) -> __m128d {
    unsafe { sse2_max_sd(a, b) }
}
/// mfence
#[inline]
pub fn _mm_mfence() {
    unsafe { sse2_mfence() }
}
/// pminsw
#[inline]
pub fn _mm_min_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// pminub
#[inline]
pub fn _mm_min_epu8(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// minpd
#[inline]
pub fn _mm_min_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// minsd
#[inline]
pub fn _mm_min_sd(a: __m128d, b: __m128d) -> __m128d {
    unsafe { sse2_min_sd(a, b) }
}
/// movq
#[inline]
pub fn _mm_move_epi64(a: __m128i) -> __m128i {
    __m128i::new(a.extract(0), 0)
}
/// movsd
#[inline]
pub fn _mm_move_sd(a: __m128d, b: __m128d) -> __m128d {
    a.replace(0, b.extract(0))
}
/// pmovmskb
#[inline]
pub fn _mm_movemask_epi8(a: __m128i) -> i32 {
    a.as_u8x16().move_mask() as i32
}
/// movmskpd
#[inline]
pub fn _mm_movemask_pd(a: __m128d) -> i32 {
    a.move_mask() as i32
}
/// pmuludq
#[inline]
pub fn _mm_mul_epu32(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// mulpd
#[inline]
pub fn _mm_mul_pd(a: __m128d, b: __m128d) -> __m128d {
    a * b
}
/// mulsd
#[inline]
pub fn _mm_mul_sd(a: __m128d, b: __m128d) -> __m128d {
    a.replace(0, a.extract(0) * b.extract(0))
}
/// pmulhw
#[inline]
pub fn _mm_mulhi_epi16(a: __m128i, b: __m128i) -> __m128i {
    a.as_i16x8().mulhi(b.as_i16x8()).as_i64x2()
}
/// pmulhuw
#[inline]
pub fn _mm_mulhi_epu16(a: __m128i, b: __m128i) -> __m128i {
    a.as_u16x8().mulhi(b.as_u16x8()).as_i64x2()
}
/// pmullw
#[inline]
pub fn _mm_mullo_epi16(a: __m128i, b: __m128i) -> __m128i {
    (a.as_i16x8() * b.as_i16x8()).as_i64x2()
}
/// orpd
#[inline]
pub fn _mm_or_pd(a: __m128d, b: __m128d) -> __m128d {
    (a.as_i64x2() | b.as_i64x2()).as_f64x2()
}
/// por
#[inline]
pub fn _mm_or_si128(a: __m128i, b: __m128i) -> __m128i {
    a | b
}
/// packsswb
#[inline]
pub fn _mm_packs_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// packssdw
#[inline]
pub fn _mm_packs_epi32(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// packuswb
#[inline]
pub fn _mm_packus_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// pause
#[inline]
pub fn _mm_pause() {
    unimplemented!()
}
/// psadbw
#[inline]
pub fn _mm_sad_epu8(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
#[inline]
pub fn _mm_set_epi16(e7: i16, e6: i16, e5: i16, e4: i16, e3: i16, e2: i16, e1: i16, e0: i16) -> __m128i {
    unimplemented!()
}
#[inline]
pub fn _mm_set_epi32(e3: i32, e2: i32, e1: i32, e0: i32) -> __m128i {
    unimplemented!()
}
#[inline]
pub fn _mm_set_epi64x(e1: i64, e0: i64) -> __m128i {
    unimplemented!()
}
#[inline]
pub fn _mm_set_epi8(e15: i8, e14: i8, e13: i8, e12: i8, e11: i8, e10: i8, e9: i8, e8: i8, e7: i8, e6: i8, e5: i8, e4: i8, e3: i8, e2: i8, e1: i8, e0: i8) -> __m128i {
    unimplemented!()
}
#[inline]
pub fn _mm_set_pd(e1: f64, e0: f64) -> __m128d {
    __m128d::new(e0, e1)
}
#[inline]
pub fn _mm_set_pd1(a: f64) -> __m128d {
    __m128d::new(a, a)
}
#[inline]
pub fn _mm_set_sd(a: f64) -> __m128d {
    __m128d::new(a, 0.)
}
#[inline]
pub fn _mm_set1_epi16(a: i16) -> __m128i {
    unimplemented!()
}
#[inline]
pub fn _mm_set1_epi32(a: i32) -> __m128i {
    unimplemented!()
}
#[inline]
pub fn _mm_set1_epi64x(a: i64) -> __m128i {
    unimplemented!()
}
#[inline]
pub fn _mm_set1_epi8(a: i8) -> __m128i {
    unimplemented!()
}
#[inline]
pub fn _mm_set1_pd(a: f64) -> __m128d {
    unimplemented!()
}
#[inline]
pub fn _mm_setr_epi16(e7: i16, e6: i16, e5: i16, e4: i16, e3: i16, e2: i16, e1: i16, e0: i16) -> __m128i {
    unimplemented!()
}
#[inline]
pub fn _mm_setr_epi32(e3: i32, e2: i32, e1: i32, e0: i32) -> __m128i {
    unimplemented!()
}
#[inline]
pub fn _mm_setr_epi8(e15: i8, e14: i8, e13: i8, e12: i8, e11: i8, e10: i8, e9: i8, e8: i8, e7: i8, e6: i8, e5: i8, e4: i8, e3: i8, e2: i8, e1: i8, e0: i8) -> __m128i {
    unimplemented!()
}
#[inline]
pub fn _mm_setr_pd(e1: f64, e0: f64) -> __m128d {
    unimplemented!()
}
/// xorpd
#[inline]
pub fn _mm_setzero_pd() -> __m128d {
    __m128d::splat(0.)
}
/// pxor
#[inline]
pub fn _mm_setzero_si128() -> __m128i {
    __m128i::splat(0)
}
/// pshufd
#[inline]
pub fn _mm_shuffle_epi32(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
}
/// shufpd
#[inline]
pub fn _mm_shuffle_pd(a: __m128d, b: __m128d, imm8: i32) -> __m128d {
    unimplemented!()
}
/// pshufhw
#[inline]
pub fn _mm_shufflehi_epi16(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
}
/// pshuflw
#[inline]
pub fn _mm_shufflelo_epi16(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
}
/// psllw
#[inline]
pub fn _mm_sll_epi16(a: __m128i, count: __m128i) -> __m128i {
    unsafe { sse2_psll_w(a.as_i16x8(), count.as_i16x8()).as_i64x2() }
}
/// pslld
#[inline]
pub fn _mm_sll_epi32(a: __m128i, count: __m128i) -> __m128i {
    unsafe { sse2_psll_d(a.as_i32x4(), count.as_i32x4()).as_i64x2() }
}
/// psllq
#[inline]
pub fn _mm_sll_epi64(a: __m128i, count: __m128i) -> __m128i {
    unsafe { sse2_psll_q(a, count) }
}
/// psllw
#[inline]
pub fn _mm_slli_epi16(a: __m128i, imm8: i32) -> __m128i {
    unsafe { sse2_pslli_w(a.as_i16x8(), imm8).as_i64x2() }
}
/// pslld
#[inline]
pub fn _mm_slli_epi32(a: __m128i, imm8: i32) -> __m128i {
    unsafe { sse2_pslli_d(a.as_i32x4(), imm8).as_i64x2() }
}
/// psllq
#[inline]
pub fn _mm_slli_epi64(a: __m128i, imm8: i32) -> __m128i {
    unsafe { sse2_pslli_q(a, imm8) }
}
/// pslldq
#[inline]
pub fn _mm_slli_si128(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
}
/// sqrtpd
#[inline]
pub fn _mm_sqrt_pd(a: __m128d) -> __m128d {
    a.sqrt()
}
/// sqrtsd
#[inline]
pub fn _mm_sqrt_sd(a: __m128d, b: __m128d) -> __m128d {
    unsafe { sse2_sqrt_sd(a) }
}
/// psraw
#[inline]
pub fn _mm_sra_epi16(a: __m128i, count: __m128i) -> __m128i {
    unsafe { sse2_psra_w(a.as_i16x8(), count.as_i16x8()).as_i64x2() }
}
/// psrad
#[inline]
pub fn _mm_sra_epi32(a: __m128i, count: __m128i) -> __m128i {
    unsafe { sse2_psra_d(a.as_i32x4(), count.as_i32x4()).as_i64x2() }
}
/// psraw
#[inline]
pub fn _mm_srai_epi16(a: __m128i, imm8: i32) -> __m128i {
    unsafe { sse2_psrai_w(a.as_i16x8(), imm8).as_i64x2() }
}
/// psrad
#[inline]
pub fn _mm_srai_epi32(a: __m128i, imm8: i32) -> __m128i {
    unsafe { sse2_psrai_d(a.as_i32x4(), imm8).as_i64x2() }
}
/// psrlw
#[inline]
pub fn _mm_srl_epi16(a: __m128i, count: __m128i) -> __m128i {
    unsafe { sse2_psrl_w(a.as_i16x8(), count.as_i16x8()).as_i64x2() }
}
/// psrld
#[inline]
pub fn _mm_srl_epi32(a: __m128i, count: __m128i) -> __m128i {
    unsafe { sse2_psrl_d(a.as_i32x4(), count.as_i32x4()).as_i64x2() }
}
/// psrlq
#[inline]
pub fn _mm_srl_epi64(a: __m128i, count: __m128i) -> __m128i {
    unsafe { sse2_psrl_q(a, count) }
}
/// psrlw
#[inline]
pub fn _mm_srli_epi16(a: __m128i, imm8: i32) -> __m128i {
    unsafe { sse2_psrli_w(a.as_i16x8(), imm8).as_i64x2() }
}
/// psrld
#[inline]
pub fn _mm_srli_epi32(a: __m128i, imm8: i32) -> __m128i {
    unsafe { sse2_psrli_d(a.as_i32x4(), imm8).as_i64x2() }
}
/// psrlq
#[inline]
pub fn _mm_srli_epi64(a: __m128i, imm8: i32) -> __m128i {
    unsafe { sse2_psrli_q(a, imm8) }
}
/// psrldq
#[inline]
pub fn _mm_srli_si128(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
}
/// movapd
#[inline]
pub unsafe fn _mm_store_pd(mem_addr: *mut f64, a: __m128d) {
    let mem_addr = mem_addr as *mut __m128d;
    *mem_addr = a;
}
#[inline]
pub unsafe fn _mm_store_pd1(mem_addr: *mut f64, a: __m128d) {
    _mm_store1_pd(mem_addr, a)
}
/// movsd
#[inline]
pub unsafe fn _mm_store_sd(mem_addr: *mut f64, a: __m128d) {
    // Cast to u8 to avoid assumptions about alignment.
    let mem_addr = mem_addr as *mut u8;
    let pa = &a as *const __m128d as *const u8;
    copy_nonoverlapping(pa, mem_addr, 8)
}
/// movdqa
#[inline]
pub unsafe fn _mm_store_si128(mem_addr: *mut __m128i, a: __m128i) {
    *mem_addr = a;
}
#[inline]
pub unsafe fn _mm_store1_pd(mem_addr: *mut f64, a: __m128d) {
    _mm_store_pd(mem_addr, __m128d::new(a.extract(0), a.extract(0)))
}
/// movhpd
#[inline]
pub unsafe fn _mm_storeh_pd(mem_addr: *mut f64, a: __m128d) {
    // Cast to u8 to avoid assumptions about alignment.
    let mem_addr = mem_addr as *mut u8;
    let pa = &a as *const __m128d as *const u8;
    copy_nonoverlapping(pa.offset(0), mem_addr, 8)
}
/// movq
#[inline]
pub unsafe fn _mm_storel_epi64(mem_addr: *mut __m128i, a: __m128i) {
    // Cast to u8 to avoid assumptions about alignment.
    let mem_addr = mem_addr as *mut u8;
    let pa = &a as *const __m128i as *const u8;
    copy_nonoverlapping(pa, mem_addr, 8)
}
/// movlpd
#[inline]
pub unsafe fn _mm_storel_pd(mem_addr: *mut f64, a: __m128d) {
    // Cast to u8 to avoid assumptions about alignment.
    let mem_addr = mem_addr as *mut u8;
    let pa = &a as *const __m128d as *const u8;
    copy_nonoverlapping(pa, mem_addr, 8)
}
#[inline]
pub unsafe fn _mm_storer_pd(mem_addr: *mut f64, a: __m128d) {
    _mm_store_pd(mem_addr, __m128d::new(a.extract(1), a.extract(0)));
}
/// movupd
#[inline]
pub unsafe fn _mm_storeu_pd(mem_addr: *mut f64, a: __m128d) {
    // Cast to u8 to avoid assumptions about alignment.
    let mem_addr = mem_addr as *mut u8;
    let pa = &a as *const __m128d as *const u8;
    copy_nonoverlapping(pa, mem_addr, 16)
}
/// movdqu
#[inline]
pub unsafe fn _mm_storeu_si128(mem_addr: *mut __m128i, a: __m128i) {
    // Cast to u8 to avoid assumptions about alignment.
    let mem_addr = mem_addr as *mut u8;
    let pa = &a as *const __m128i as *const u8;
    copy_nonoverlapping(pa, mem_addr, 16)
}
/// psubw
#[inline]
pub fn _mm_sub_epi16(a: __m128i, b: __m128i) -> __m128i {
    (a.as_i16x8() + b.as_i16x8()).as_i64x2()
}
/// psubd
#[inline]
pub fn _mm_sub_epi32(a: __m128i, b: __m128i) -> __m128i {
    (a.as_i32x4() - b.as_i32x4()).as_i64x2()
}
/// psubq
#[inline]
pub fn _mm_sub_epi64(a: __m128i, b: __m128i) -> __m128i {
    a - b
}
/// psubb
#[inline]
pub fn _mm_sub_epi8(a: __m128i, b: __m128i) -> __m128i {
    (a.as_i8x16() + b.as_i8x16()).as_i64x2()
}
/// subpd
#[inline]
pub fn _mm_sub_pd(a: __m128d, b: __m128d) -> __m128d {
    a - b
}
/// subsd
#[inline]
pub fn _mm_sub_sd(a: __m128d, b: __m128d) -> __m128d {
    a.replace(0, a.extract(0) - b.extract(0))
}
/// psubsw
#[inline]
pub fn _mm_subs_epi16(a: __m128i, b: __m128i) -> __m128i {
    a.as_i16x8().subs(b.as_i16x8()).as_i64x2()
}
/// psubsb
#[inline]
pub fn _mm_subs_epi8(a: __m128i, b: __m128i) -> __m128i {
    a.as_i8x16().subs(b.as_i8x16()).as_i64x2()
}
/// psubusw
#[inline]
pub fn _mm_subs_epu16(a: __m128i, b: __m128i) -> __m128i {
    a.as_u16x8().subs(b.as_u16x8()).as_i64x2()
}
/// psubusb
#[inline]
pub fn _mm_subs_epu8(a: __m128i, b: __m128i) -> __m128i {
    a.as_u8x16().subs(b.as_u8x16()).as_i64x2()
}
/// ucomisd
#[inline]
pub fn _mm_ucomieq_sd(a: __m128d, b: __m128d) -> i32 {
    unsafe { sse2_ucomieq_sd(a, b) }
}
/// ucomisd
#[inline]
pub fn _mm_ucomige_sd(a: __m128d, b: __m128d) -> i32 {
    unsafe { sse2_ucomige_sd(a, b) }
}
/// ucomisd
#[inline]
pub fn _mm_ucomigt_sd(a: __m128d, b: __m128d) -> i32 {
    unsafe { sse2_ucomigt_sd(a, b) }
}
/// ucomisd
#[inline]
pub fn _mm_ucomile_sd(a: __m128d, b: __m128d) -> i32 {
    unsafe { sse2_ucomile_sd(a, b) }
}
/// ucomisd
#[inline]
pub fn _mm_ucomilt_sd(a: __m128d, b: __m128d) -> i32 {
    unsafe { sse2_ucomilt_sd(a, b) }
}
/// ucomisd
#[inline]
pub fn _mm_ucomineq_sd(a: __m128d, b: __m128d) -> i32 {
    unsafe { sse2_ucomineq_sd(a, b) }
}
/// punpckhwd
#[inline]
pub fn _mm_unpackhi_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// punpckhdq
#[inline]
pub fn _mm_unpackhi_epi32(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// punpckhqdq
#[inline]
pub fn _mm_unpackhi_epi64(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// punpckhbw
#[inline]
pub fn _mm_unpackhi_epi8(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// unpckhpd
#[inline]
pub fn _mm_unpackhi_pd(a: __m128d, b: __m128d) -> __m128d {
    __m128d::new(a.extract(1), b.extract(1))
}
/// punpcklwd
#[inline]
pub fn _mm_unpacklo_epi16(a: __m128i, b: __m128i) -> __m128i {
    let a = a.as_i16x8();
    let b = b.as_i16x8();
    let r: i16x8 = unsafe { simd_shuffle8(a, b, [0, 8, 1, 9, 2, 10, 3, 11]) };
    r.as_i64x2()
}
/// punpckldq
#[inline]
pub fn _mm_unpacklo_epi32(a: __m128i, b: __m128i) -> __m128i {
    let r: i32x4 = unsafe { simd_shuffle4(a.as_i32x4(), b.as_i32x4(), [0, 4, 1, 5]) };
    r.as_i64x2()
}
/// punpcklqdq
#[inline]
pub fn _mm_unpacklo_epi64(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// punpcklbw
#[inline]
pub fn _mm_unpacklo_epi8(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// unpcklpd
#[inline]
pub fn _mm_unpacklo_pd(a: __m128d, b: __m128d) -> __m128d {
    __m128d::new(a.extract(0), b.extract(0))
}
/// xorpd
#[inline]
pub fn _mm_xor_pd(a: __m128d, b: __m128d) -> __m128d {
    (a.as_i64x2() ^ b.as_i64x2()).as_f64x2()
}
/// pxor
#[inline]
pub fn _mm_xor_si128(a: __m128i, b: __m128i) -> __m128i {
    a ^ b
}

/// The methods in this module can't be implemented because Rust doesn't
/// expose nontemporal loads.
pub mod unimplemented_nontemporal {
    use __m128d;
    use __m128i;

    /// movntpd
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_stream_pd(mem_addr: *mut f64, a: __m128d) {
        unimplemented!()
    }
    /// movntdq
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_stream_si128(mem_addr: *mut __m128i, a: __m128i) {
        unimplemented!()
    }
    /// movnti
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_stream_si32(mem_addr: *mut i32, a: i32) {
        unimplemented!()
    }
    /// movnti
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_stream_si64(mem_addr: *mut i64, a: i64) {
        unimplemented!()
    }
}

/// The methods in this module can't be implemented because Rust doesn't
/// expose the LLVM x86_mmx type.
pub mod unimplemented_mmx {
    use __m64;
    use __m128i;
    use __m128d;

    /// paddq
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_add_si64(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// cvtpd2pi
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_cvtpd_pi32(a: __m128d) -> __m64 {
        unimplemented!()
    }
    /// cvtpi2pd
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_cvtpi32_pd(a: __m64) -> __m128d {
        unimplemented!()
    }
    /// cvttpd2pi
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_cvttpd_pi32(a: __m128d) -> __m64 {
        unimplemented!()
    }
    /// movdq2q
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_movepi64_pi64(a: __m128i) -> __m64 {
        unimplemented!()
    }
    /// movq2dq
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_movpi64_epi64(a: __m64) -> __m128i {
        unimplemented!()
    }
    /// pmuludq
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_mul_su32(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// Not yet implemented.
    #[inline]
    pub fn _mm_set_epi64(e1: __m64, e0: __m64) -> __m128i {
        unimplemented!()
    }
    /// Not yet implemented.
    #[inline]
    pub fn _mm_set1_epi64(a: __m64) -> __m128i {
        unimplemented!()
    }
    /// Not yet implemented.
    #[inline]
    pub fn _mm_setr_epi64(e1: __m64, e0: __m64) -> __m128i {
        unimplemented!()
    }
    /// psubq
    ///
    /// Not yet implemented.
    #[inline]
    pub fn _mm_sub_si64(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
}
