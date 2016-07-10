// Copyright 2016 Eli Friedman.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#![allow(unused_variables)]

use conversions::Convert128;
use __m64;
use __m128;
use __m128i;
use __m128d;

/// paddw
#[inline]
pub fn _mm_add_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// paddd
#[inline]
pub fn _mm_add_epi32(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// paddq
#[inline]
pub fn _mm_add_epi64(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// paddb
#[inline]
pub fn _mm_add_epi8(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// addpd
#[inline]
pub fn _mm_add_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// addsd
#[inline]
pub fn _mm_add_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// paddq
#[inline]
pub fn _mm_add_si64(a: __m64, b: __m64) -> __m64 {
    unimplemented!()
}
/// paddsw
#[inline]
pub fn _mm_adds_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// paddsb
#[inline]
pub fn _mm_adds_epi8(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// paddusw
#[inline]
pub fn _mm_adds_epu16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// paddusb
#[inline]
pub fn _mm_adds_epu8(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// andpd
#[inline]
pub fn _mm_and_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// pand
#[inline]
pub fn _mm_and_si128(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// andnpd
#[inline]
pub fn _mm_andnot_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// pandn
#[inline]
pub fn _mm_andnot_si128(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// pavgw
#[inline]
pub fn _mm_avg_epu16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// pavgb
#[inline]
pub fn _mm_avg_epu8(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// pslldq
#[inline]
pub fn _mm_bslli_si128(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
}
/// psrldq
#[inline]
pub fn _mm_bsrli_si128(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
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
    unimplemented!()
}
/// pcmpeqw
#[inline]
pub fn _mm_cmpeq_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// pcmpeqd
#[inline]
pub fn _mm_cmpeq_epi32(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// pcmpeqb
#[inline]
pub fn _mm_cmpeq_epi8(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpeq_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmpsd
#[inline]
pub fn _mm_cmpeq_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpge_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmpsd
#[inline]
pub fn _mm_cmpge_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// pcmpgtw
#[inline]
pub fn _mm_cmpgt_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// pcmpgtd
#[inline]
pub fn _mm_cmpgt_epi32(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// pcmpgtb
#[inline]
pub fn _mm_cmpgt_epi8(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpgt_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmpsd
#[inline]
pub fn _mm_cmpgt_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmple_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmpsd
#[inline]
pub fn _mm_cmple_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// pcmpgtw
#[inline]
pub fn _mm_cmplt_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// pcmpgtd
#[inline]
pub fn _mm_cmplt_epi32(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// pcmpgtb
#[inline]
pub fn _mm_cmplt_epi8(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmplt_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmpsd
#[inline]
pub fn _mm_cmplt_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpneq_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmpsd
#[inline]
pub fn _mm_cmpneq_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpnge_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmpsd
#[inline]
pub fn _mm_cmpnge_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpngt_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmpsd
#[inline]
pub fn _mm_cmpngt_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpnle_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmpsd
#[inline]
pub fn _mm_cmpnle_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpnlt_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmpsd
#[inline]
pub fn _mm_cmpnlt_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpord_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmpsd
#[inline]
pub fn _mm_cmpord_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmppd
#[inline]
pub fn _mm_cmpunord_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// cmpsd
#[inline]
pub fn _mm_cmpunord_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// comisd
#[inline]
pub fn _mm_comieq_sd(a: __m128d, b: __m128d) -> i32 {
    unimplemented!()
}
/// comisd
#[inline]
pub fn _mm_comige_sd(a: __m128d, b: __m128d) -> i32 {
    unimplemented!()
}
/// comisd
#[inline]
pub fn _mm_comigt_sd(a: __m128d, b: __m128d) -> i32 {
    unimplemented!()
}
/// comisd
#[inline]
pub fn _mm_comile_sd(a: __m128d, b: __m128d) -> i32 {
    unimplemented!()
}
/// comisd
#[inline]
pub fn _mm_comilt_sd(a: __m128d, b: __m128d) -> i32 {
    unimplemented!()
}
/// comisd
#[inline]
pub fn _mm_comineq_sd(a: __m128d, b: __m128d) -> i32 {
    unimplemented!()
}
/// cvtdq2pd
#[inline]
pub fn _mm_cvtepi32_pd(a: __m128i) -> __m128d {
    unimplemented!()
}
/// cvtdq2ps
#[inline]
pub fn _mm_cvtepi32_ps(a: __m128i) -> __m128 {
    unimplemented!()
}
/// cvtpd2dq
#[inline]
pub fn _mm_cvtpd_epi32(a: __m128d) -> __m128i {
    unimplemented!()
}
/// cvtpd2pi
#[inline]
pub fn _mm_cvtpd_pi32(a: __m128d) -> __m64 {
    unimplemented!()
}
/// cvtpd2ps
#[inline]
pub fn _mm_cvtpd_ps(a: __m128d) -> __m128 {
    unimplemented!()
}
/// cvtpi2pd
#[inline]
pub fn _mm_cvtpi32_pd(a: __m64) -> __m128d {
    unimplemented!()
}
/// cvtps2dq
#[inline]
pub fn _mm_cvtps_epi32(a: __m128) -> __m128i {
    unimplemented!()
}
/// cvtps2pd
#[inline]
pub fn _mm_cvtps_pd(a: __m128) -> __m128d {
    unimplemented!()
}
/// movsd
#[inline]
pub fn _mm_cvtsd_f64(a: __m128d) -> f64 {
    unimplemented!()
}
/// cvtsd2si
#[inline]
pub fn _mm_cvtsd_si32(a: __m128d) -> i32 {
    unimplemented!()
}
/// cvtsd2si
#[inline]
pub fn _mm_cvtsd_si64(a: __m128d) -> i64 {
    unimplemented!()
}
/// cvtsd2si
#[inline]
pub fn _mm_cvtsd_si64x(a: __m128d) -> i64 {
    unimplemented!()
}
/// cvtsd2ss
#[inline]
pub fn _mm_cvtsd_ss(a: __m128, b: __m128d) -> __m128 {
    unimplemented!()
}
/// movd
#[inline]
pub fn _mm_cvtsi128_si32(a: __m128i) -> i32 {
    unimplemented!()
}
/// movq
#[inline]
pub fn _mm_cvtsi128_si64(a: __m128i) -> i64 {
    unimplemented!()
}
/// movq
#[inline]
pub fn _mm_cvtsi128_si64x(a: __m128i) -> i64 {
    unimplemented!()
}
/// cvtsi2sd
#[inline]
pub fn _mm_cvtsi32_sd(a: __m128d, b: i32) -> __m128d {
    unimplemented!()
}
/// movd
#[inline]
pub fn _mm_cvtsi32_si128(a: i32) -> __m128i {
    unimplemented!()
}
/// cvtsi2sd
#[inline]
pub fn _mm_cvtsi64_sd(a: __m128d, b: i64) -> __m128d {
    unimplemented!()
}
/// movq
#[inline]
pub fn _mm_cvtsi64_si128(a: i64) -> __m128i {
    unimplemented!()
}
/// cvtsi2sd
#[inline]
pub fn _mm_cvtsi64x_sd(a: __m128d, b: i64) -> __m128d {
    unimplemented!()
}
/// movq
#[inline]
pub fn _mm_cvtsi64x_si128(a: i64) -> __m128i {
    unimplemented!()
}
/// cvtss2sd
#[inline]
pub fn _mm_cvtss_sd(a: __m128d, b: __m128) -> __m128d {
    unimplemented!()
}
/// cvttpd2dq
#[inline]
pub fn _mm_cvttpd_epi32(a: __m128d) -> __m128i {
    unimplemented!()
}
/// cvttpd2pi
#[inline]
pub fn _mm_cvttpd_pi32(a: __m128d) -> __m64 {
    unimplemented!()
}
/// cvttps2dq
#[inline]
pub fn _mm_cvttps_epi32(a: __m128) -> __m128i {
    unimplemented!()
}
/// cvttsd2si
#[inline]
pub fn _mm_cvttsd_si32(a: __m128d) -> i32 {
    unimplemented!()
}
/// cvttsd2si
#[inline]
pub fn _mm_cvttsd_si64(a: __m128d) -> i64 {
    unimplemented!()
}
/// cvttsd2si
#[inline]
pub fn _mm_cvttsd_si64x(a: __m128d) -> i64 {
    unimplemented!()
}
/// divpd
#[inline]
pub fn _mm_div_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// divsd
#[inline]
pub fn _mm_div_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// pextrw
#[inline]
pub fn _mm_extract_epi16(a: __m128i, imm8: i32) -> i32 {
    unimplemented!()
}
/// pinsrw
#[inline]
pub fn _mm_insert_epi16(a: __m128i, i: i32, imm8: i32) -> __m128i {
    unimplemented!()
}
/// lfence
#[inline]
pub fn _mm_lfence() {
    unimplemented!()
}
/// movapd
#[inline]
pub fn _mm_load_pd(mem_addr: *const f64) -> __m128d {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_load_pd1(mem_addr: *const f64) -> __m128d {
    unimplemented!()
}
/// movsd
#[inline]
pub fn _mm_load_sd(mem_addr: *const f64) -> __m128d {
    unimplemented!()
}
/// movdqa
#[inline]
pub fn _mm_load_si128(mem_addr: *const __m128i) -> __m128i {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_load1_pd(mem_addr: *const f64) -> __m128d {
    unimplemented!()
}
/// movhpd
#[inline]
pub fn _mm_loadh_pd(a: __m128d, mem_addr: *const f64) -> __m128d {
    unimplemented!()
}
/// movq
#[inline]
pub fn _mm_loadl_epi64(mem_addr: *const __m128i) -> __m128i {
    unimplemented!()
}
/// movlpd
#[inline]
pub fn _mm_loadl_pd(a: __m128d, mem_addr: *const f64) -> __m128d {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_loadr_pd(mem_addr: *const f64) -> __m128d {
    unimplemented!()
}
/// movupd
#[inline]
pub fn _mm_loadu_pd(mem_addr: *const f64) -> __m128d {
    unimplemented!()
}
/// movdqu
#[inline]
pub fn _mm_loadu_si128(mem_addr: *const __m128i) -> __m128i {
    unimplemented!()
}
/// pmaddwd
#[inline]
pub fn _mm_madd_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// maskmovdqu
#[inline]
pub fn _mm_maskmoveu_si128(a: __m128i, mask: __m128i, mem_addr: *mut i8) {
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
    unimplemented!()
}
/// maxsd
#[inline]
pub fn _mm_max_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// mfence
#[inline]
pub fn _mm_mfence() {
    unimplemented!()
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
    unimplemented!()
}
/// movq
#[inline]
pub fn _mm_move_epi64(a: __m128i) -> __m128i {
    unimplemented!()
}
/// movsd
#[inline]
pub fn _mm_move_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// pmovmskb
#[inline]
pub fn _mm_movemask_epi8(a: __m128i) -> i32 {
    unimplemented!()
}
/// movmskpd
#[inline]
pub fn _mm_movemask_pd(a: __m128d) -> i32 {
    unimplemented!()
}
/// movdq2q
#[inline]
pub fn _mm_movepi64_pi64(a: __m128i) -> __m64 {
    unimplemented!()
}
/// movq2dq
#[inline]
pub fn _mm_movpi64_epi64(a: __m64) -> __m128i {
    unimplemented!()
}
/// pmuludq
#[inline]
pub fn _mm_mul_epu32(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// mulpd
#[inline]
pub fn _mm_mul_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// mulsd
#[inline]
pub fn _mm_mul_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// pmuludq
#[inline]
pub fn _mm_mul_su32(a: __m64, b: __m64) -> __m64 {
    unimplemented!()
}
/// pmulhw
#[inline]
pub fn _mm_mulhi_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// pmulhuw
#[inline]
pub fn _mm_mulhi_epu16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// pmullw
#[inline]
pub fn _mm_mullo_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// orpd
#[inline]
pub fn _mm_or_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// por
#[inline]
pub fn _mm_or_si128(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
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
/// ...
#[inline]
pub fn _mm_set_epi16(e7: i16, e6: i16, e5: i16, e4: i16, e3: i16, e2: i16, e1: i16, e0: i16) -> __m128i {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_set_epi32(e3: i32, e2: i32, e1: i32, e0: i32) -> __m128i {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_set_epi64(e1: __m64, e0: __m64) -> __m128i {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_set_epi64x(e1: i64, e0: i64) -> __m128i {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_set_epi8(e15: i8, e14: i8, e13: i8, e12: i8, e11: i8, e10: i8, e9: i8, e8: i8, e7: i8, e6: i8, e5: i8, e4: i8, e3: i8, e2: i8, e1: i8, e0: i8) -> __m128i {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_set_pd(e1: f64, e0: f64) -> __m128d {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_set_pd1(a: f64) -> __m128d {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_set_sd(a: f64) -> __m128d {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_set1_epi16(a: i16) -> __m128i {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_set1_epi32(a: i32) -> __m128i {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_set1_epi64(a: __m64) -> __m128i {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_set1_epi64x(a: i64) -> __m128i {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_set1_epi8(a: i8) -> __m128i {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_set1_pd(a: f64) -> __m128d {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_setr_epi16(e7: i16, e6: i16, e5: i16, e4: i16, e3: i16, e2: i16, e1: i16, e0: i16) -> __m128i {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_setr_epi32(e3: i32, e2: i32, e1: i32, e0: i32) -> __m128i {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_setr_epi64(e1: __m64, e0: __m64) -> __m128i {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_setr_epi8(e15: i8, e14: i8, e13: i8, e12: i8, e11: i8, e10: i8, e9: i8, e8: i8, e7: i8, e6: i8, e5: i8, e4: i8, e3: i8, e2: i8, e1: i8, e0: i8) -> __m128i {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_setr_pd(e1: f64, e0: f64) -> __m128d {
    unimplemented!()
}
/// xorpd
#[inline]
pub fn _mm_setzero_pd() -> __m128d {
    unimplemented!()
}
/// pxor
#[inline]
pub fn _mm_setzero_si128() -> __m128i {
    unimplemented!()
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
    unimplemented!()
}
/// pslld
#[inline]
pub fn _mm_sll_epi32(a: __m128i, count: __m128i) -> __m128i {
    unimplemented!()
}
/// psllq
#[inline]
pub fn _mm_sll_epi64(a: __m128i, count: __m128i) -> __m128i {
    unimplemented!()
}
/// psllw
#[inline]
pub fn _mm_slli_epi16(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
}
/// pslld
#[inline]
pub fn _mm_slli_epi32(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
}
/// psllq
#[inline]
pub fn _mm_slli_epi64(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
}
/// pslldq
#[inline]
pub fn _mm_slli_si128(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
}
/// sqrtpd
#[inline]
pub fn _mm_sqrt_pd(a: __m128d) -> __m128d {
    unimplemented!()
}
/// sqrtsd
#[inline]
pub fn _mm_sqrt_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// psraw
#[inline]
pub fn _mm_sra_epi16(a: __m128i, count: __m128i) -> __m128i {
    unimplemented!()
}
/// psrad
#[inline]
pub fn _mm_sra_epi32(a: __m128i, count: __m128i) -> __m128i {
    unimplemented!()
}
/// psraw
#[inline]
pub fn _mm_srai_epi16(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
}
/// psrad
#[inline]
pub fn _mm_srai_epi32(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
}
/// psrlw
#[inline]
pub fn _mm_srl_epi16(a: __m128i, count: __m128i) -> __m128i {
    unimplemented!()
}
/// psrld
#[inline]
pub fn _mm_srl_epi32(a: __m128i, count: __m128i) -> __m128i {
    unimplemented!()
}
/// psrlq
#[inline]
pub fn _mm_srl_epi64(a: __m128i, count: __m128i) -> __m128i {
    unimplemented!()
}
/// psrlw
#[inline]
pub fn _mm_srli_epi16(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
}
/// psrld
#[inline]
pub fn _mm_srli_epi32(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
}
/// psrlq
#[inline]
pub fn _mm_srli_epi64(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
}
/// psrldq
#[inline]
pub fn _mm_srli_si128(a: __m128i, imm8: i32) -> __m128i {
    unimplemented!()
}
/// movapd
#[inline]
pub fn _mm_store_pd(mem_addr: *mut f64, a: __m128d) {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_store_pd1(mem_addr: *mut f64, a: __m128d) {
    unimplemented!()
}
/// movsd
#[inline]
pub fn _mm_store_sd(mem_addr: *mut f64, a: __m128d) {
    unimplemented!()
}
/// movdqa
#[inline]
pub fn _mm_store_si128(mem_addr: *mut __m128i, a: __m128i) {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_store1_pd(mem_addr: *mut f64, a: __m128d) {
    unimplemented!()
}
/// movhpd
#[inline]
pub fn _mm_storeh_pd(mem_addr: *mut f64, a: __m128d) {
    unimplemented!()
}
/// movq
#[inline]
pub fn _mm_storel_epi64(mem_addr: *mut __m128i, a: __m128i) {
    unimplemented!()
}
/// movlpd
#[inline]
pub fn _mm_storel_pd(mem_addr: *mut f64, a: __m128d) {
    unimplemented!()
}
/// ...
#[inline]
pub fn _mm_storer_pd(mem_addr: *mut f64, a: __m128d) {
    unimplemented!()
}
/// movupd
#[inline]
pub fn _mm_storeu_pd(mem_addr: *mut f64, a: __m128d) {
    unimplemented!()
}
/// movdqu
#[inline]
pub fn _mm_storeu_si128(mem_addr: *mut __m128i, a: __m128i) {
    unimplemented!()
}
/// movntpd
#[inline]
pub fn _mm_stream_pd(mem_addr: *mut f64, a: __m128d) {
    unimplemented!()
}
/// movntdq
#[inline]
pub fn _mm_stream_si128(mem_addr: *mut __m128i, a: __m128i) {
    unimplemented!()
}
/// movnti
#[inline]
pub fn _mm_stream_si32(mem_addr: *mut i32, a: i32) {
    unimplemented!()
}
/// movnti
#[inline]
pub fn _mm_stream_si64(mem_addr: *mut i64, a: i64) {
    unimplemented!()
}
/// psubw
#[inline]
pub fn _mm_sub_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// psubd
#[inline]
pub fn _mm_sub_epi32(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// psubq
#[inline]
pub fn _mm_sub_epi64(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// psubb
#[inline]
pub fn _mm_sub_epi8(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// subpd
#[inline]
pub fn _mm_sub_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// subsd
#[inline]
pub fn _mm_sub_sd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// psubq
#[inline]
pub fn _mm_sub_si64(a: __m64, b: __m64) -> __m64 {
    unimplemented!()
}
/// psubsw
#[inline]
pub fn _mm_subs_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// psubsb
#[inline]
pub fn _mm_subs_epi8(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// psubusw
#[inline]
pub fn _mm_subs_epu16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// psubusb
#[inline]
pub fn _mm_subs_epu8(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// ucomisd
#[inline]
pub fn _mm_ucomieq_sd(a: __m128d, b: __m128d) -> i32 {
    unimplemented!()
}
/// ucomisd
#[inline]
pub fn _mm_ucomige_sd(a: __m128d, b: __m128d) -> i32 {
    unimplemented!()
}
/// ucomisd
#[inline]
pub fn _mm_ucomigt_sd(a: __m128d, b: __m128d) -> i32 {
    unimplemented!()
}
/// ucomisd
#[inline]
pub fn _mm_ucomile_sd(a: __m128d, b: __m128d) -> i32 {
    unimplemented!()
}
/// ucomisd
#[inline]
pub fn _mm_ucomilt_sd(a: __m128d, b: __m128d) -> i32 {
    unimplemented!()
}
/// ucomisd
#[inline]
pub fn _mm_ucomineq_sd(a: __m128d, b: __m128d) -> i32 {
    unimplemented!()
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
    unimplemented!()
}
/// punpcklwd
#[inline]
pub fn _mm_unpacklo_epi16(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
/// punpckldq
#[inline]
pub fn _mm_unpacklo_epi32(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
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
    unimplemented!()
}
/// xorpd
#[inline]
pub fn _mm_xor_pd(a: __m128d, b: __m128d) -> __m128d {
    unimplemented!()
}
/// pxor
#[inline]
pub fn _mm_xor_si128(a: __m128i, b: __m128i) -> __m128i {
    unimplemented!()
}
