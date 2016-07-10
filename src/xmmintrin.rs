// Copyright 2016 Eli Friedman.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#![allow(unused_variables)]

use std::mem::transmute;
use std::mem::uninitialized;
use std::ptr::copy_nonoverlapping;
use simd;
use __m128;
use __m128i;
use __m64;

fn convert_bool32fx4_to_m128(a: simd::bool32fx4) -> __m128 {
    unsafe { transmute(a) }
}
fn convert_m128i_to_m128(a: __m128i) -> __m128 {
    unsafe { transmute(a) }
}
fn convert_m128_to_m128i(a: __m128) -> __m128i {
    unsafe { transmute(a) }
}

// Declarations copied from the llvmint crate.
#[allow(improper_ctypes)]
extern {
    #[link_name = "llvm.x86.sse.sqrt.ss"]
    pub fn sse_sqrt_ss(a: __m128) -> __m128;
    #[link_name = "llvm.x86.sse.cvtss2si"]
    pub fn sse_cvtss2si(a: __m128) -> i32;
    #[link_name = "llvm.x86.sse.cvttss2si"]
    pub fn sse_cvttss2si(a: __m128) -> i32;
    #[link_name = "llvm.x86.sse.cvtss2si64"]
    pub fn sse_cvtss2si64(a: __m128) -> i64;
    #[link_name = "llvm.x86.sse.cvttss2si64"]
    pub fn sse_cvttss2si64(a: __m128) -> i64;
    #[link_name = "llvm.x86.sse.movmsk.ps"]
    pub fn sse_movmsk_ps(a: __m128) -> i32;
    #[link_name = "llvm.x86.sse.sfence"]
    pub fn sse_sfence() ->();
    #[link_name = "llvm.x86.sse.rcp.ss"]
    pub fn sse_rcp_ss(a: __m128) -> __m128;
    #[link_name = "llvm.x86.sse.rsqrt.ss"]
    pub fn sse_rsqrt_ss(a: __m128) -> __m128;
    #[link_name = "llvm.x86.sse.min.ss"]
    pub fn sse_min_ss(a: __m128, b: __m128) -> __m128;
    #[link_name = "llvm.x86.sse.max.ss"]
    pub fn sse_max_ss(a: __m128, b: __m128) -> __m128;
    #[link_name = "llvm.x86.sse.cmp.ss"]
    pub fn sse_cmp_ss(a: __m128, b: __m128, c: i8) -> __m128;
    #[link_name = "llvm.x86.sse.comieq.ss"]
    pub fn sse_comieq_ss(a: __m128, b: __m128) -> i32;
    #[link_name = "llvm.x86.sse.comilt.ss"]
    pub fn sse_comilt_ss(a: __m128, b: __m128) -> i32;
    #[link_name = "llvm.x86.sse.comile.ss"]
    pub fn sse_comile_ss(a: __m128, b: __m128) -> i32;
    #[link_name = "llvm.x86.sse.comigt.ss"]
    pub fn sse_comigt_ss(a: __m128, b: __m128) -> i32;
    #[link_name = "llvm.x86.sse.comige.ss"]
    pub fn sse_comige_ss(a: __m128, b: __m128) -> i32;
    #[link_name = "llvm.x86.sse.comineq.ss"]
    pub fn sse_comineq_ss(a: __m128, b: __m128) -> i32;
    #[link_name = "llvm.x86.sse.ucomieq.ss"]
    pub fn sse_ucomieq_ss(a: __m128, b: __m128) -> i32;
    #[link_name = "llvm.x86.sse.ucomilt.ss"]
    pub fn sse_ucomilt_ss(a: __m128, b: __m128) -> i32;
    #[link_name = "llvm.x86.sse.ucomile.ss"]
    pub fn sse_ucomile_ss(a: __m128, b: __m128) -> i32;
    #[link_name = "llvm.x86.sse.ucomigt.ss"]
    pub fn sse_ucomigt_ss(a: __m128, b: __m128) -> i32;
    #[link_name = "llvm.x86.sse.ucomige.ss"]
    pub fn sse_ucomige_ss(a: __m128, b: __m128) -> i32;
    #[link_name = "llvm.x86.sse.ucomineq.ss"]
    pub fn sse_ucomineq_ss(a: __m128, b: __m128) -> i32;
    #[link_name = "llvm.x86.sse.cvtsi2ss"]
    pub fn sse_cvtsi2ss(a: __m128, b: i32) -> __m128;
    #[link_name = "llvm.x86.sse.cvtsi642ss"]
    pub fn sse_cvtsi642ss(a: __m128, b: i64) -> __m128;
    #[link_name = "llvm.x86.sse.stmxcsr"]
    pub fn sse_stmxcsr(a: *mut i8) -> ();
    #[link_name = "llvm.x86.sse.ldmxcsr"]
    pub fn sse_ldmxcsr(a: *mut i8) -> ();
    #[link_name = "llvm.prefetch"]
    pub fn prefetch(a: *const i8, b: i32, c: i32, d: i32) -> ();
}

/// addss
#[inline]
pub fn _mm_add_ss(a: __m128, b: __m128) -> __m128 {
    a.replace(0, a.extract(0) + b.extract(0))
}
/// addps
#[inline]
pub fn _mm_add_ps(a: __m128, b: __m128) -> __m128 {
    a + b
}
/// subss
#[inline]
pub fn _mm_sub_ss(a: __m128, b: __m128) -> __m128 {
    a.replace(0, a.extract(0) - b.extract(0))
}
/// subps
#[inline]
pub fn _mm_sub_ps(a: __m128, b: __m128) -> __m128 {
    a - b
}
/// mulss
#[inline]
pub fn _mm_mul_ss(a: __m128, b: __m128) -> __m128 {
    a.replace(0, a.extract(0) * b.extract(0))
}
/// mulps
#[inline]
pub fn _mm_mul_ps(a: __m128, b: __m128) -> __m128 {
    a * b
}
/// divss
#[inline]
pub fn _mm_div_ss(a: __m128, b: __m128) -> __m128 {
    a.replace(0, a.extract(0) / b.extract(0))
}
/// divps
#[inline]
pub fn _mm_div_ps(a: __m128, b: __m128) -> __m128 {
    a / b
}
/// sqrtss
#[inline]
pub fn _mm_sqrt_ss(a: __m128) -> __m128 {
    unsafe { sse_sqrt_ss(a) }
}
/// sqrtps
#[inline]
pub fn _mm_sqrt_ps(a: __m128) -> __m128 {
    a.sqrt()
}
/// rcpss
#[inline]
pub fn _mm_rcp_ss(a: __m128) -> __m128 {
    unsafe { sse_rcp_ss(a) }
}
/// rcpps
#[inline]
pub fn _mm_rcp_ps(a: __m128) -> __m128 {
    a.approx_reciprocal()
}
/// rsqrtss
#[inline]
pub fn _mm_rsqrt_ss(a: __m128) -> __m128 {
    unsafe { sse_rsqrt_ss(a) }
}
/// rsqrtps
#[inline]
pub fn _mm_rsqrt_ps(a: __m128) -> __m128 {
    a.approx_rsqrt()
}
/// minss
#[inline]
pub fn _mm_min_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_min_ss(a, b) }
}
/// minps
#[inline]
pub fn _mm_min_ps(a: __m128, b: __m128) -> __m128 {
    a.min(b)
}
/// maxss
#[inline]
pub fn _mm_max_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_max_ss(a, b) }
}
/// maxps
#[inline]
pub fn _mm_max_ps(a: __m128, b: __m128) -> __m128 {
    a.max(b)
}
/// andps
#[inline]
pub fn _mm_and_ps(a: __m128, b: __m128) -> __m128 {
    convert_m128i_to_m128(convert_m128_to_m128i(a) & convert_m128_to_m128i(b))
}
/// andnps
#[inline]
pub fn _mm_andnot_ps(a: __m128, b: __m128) -> __m128 {
    // FIXME: Not operator from simd doesn't get inlined?!
    convert_m128i_to_m128((convert_m128_to_m128i(a) ^ __m128i::splat(!0)) & convert_m128_to_m128i(b))
}
/// orps
#[inline]
pub fn _mm_or_ps(a: __m128, b: __m128) -> __m128 {
    convert_m128i_to_m128(convert_m128_to_m128i(a) | convert_m128_to_m128i(b))
}
/// xorps
#[inline]
pub fn _mm_xor_ps(a: __m128, b: __m128) -> __m128 {
    convert_m128i_to_m128(convert_m128_to_m128i(a) ^ convert_m128_to_m128i(b))
}
/// cmpss
#[inline]
pub fn _mm_cmpeq_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_cmp_ss(a, b, 0) }
}
/// cmpps
#[inline]
pub fn _mm_cmpeq_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.eq(b))
}
/// cmpss
#[inline]
pub fn _mm_cmplt_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_cmp_ss(a, b, 1) }
}
/// cmpps
#[inline]
pub fn _mm_cmplt_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.lt(b))
}
/// cmpss
#[inline]
pub fn _mm_cmple_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_cmp_ss(a, b, 2) }
}
/// cmpps
#[inline]
pub fn _mm_cmple_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.le(b))
}
/// cmpss
#[inline]
pub fn _mm_cmpgt_ss(a: __m128, b: __m128) -> __m128 {
    _mm_move_ss(a, _mm_cmplt_ss(b, a))
}
/// cmpps
#[inline]
pub fn _mm_cmpgt_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.gt(b))
}
/// cmpss
#[inline]
pub fn _mm_cmpge_ss(a: __m128, b: __m128) -> __m128 {
    _mm_move_ss(a, _mm_cmple_ss(b, a))
}
/// cmpps
#[inline]
pub fn _mm_cmpge_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.ge(b))
}
/// cmpss
#[inline]
pub fn _mm_cmpneq_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_cmp_ss(a, b, 4) }
}
/// cmpps
#[inline]
pub fn _mm_cmpneq_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.ne(b))
}
/// cmpss
#[inline]
pub fn _mm_cmpnlt_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_cmp_ss(a, b, 5) }
}
/// cmpps
#[inline]
pub fn _mm_cmpnlt_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.ge(b) | (a.ne(a) | b.ne(b)))
}
/// cmpss
#[inline]
pub fn _mm_cmpnle_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_cmp_ss(a, b, 6) }
}
/// cmpps
#[inline]
pub fn _mm_cmpnle_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.gt(b) | (a.ne(a) | b.ne(b)))
}
/// cmpss
#[inline]
pub fn _mm_cmpngt_ss(a: __m128, b: __m128) -> __m128 {
    _mm_move_ss(a, _mm_cmpnlt_ss(b, a))
}
/// cmpps
#[inline]
pub fn _mm_cmpngt_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.le(b) | (a.ne(a) | b.ne(b)))
}
/// cmpss
#[inline]
pub fn _mm_cmpnge_ss(a: __m128, b: __m128) -> __m128 {
    _mm_move_ss(a, _mm_cmpnle_ss(b, a))
}
/// cmpps
#[inline]
pub fn _mm_cmpnge_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.lt(b) |(a.ne(a) | b.ne(b)))
}
/// cmpss
#[inline]
pub fn _mm_cmpord_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_cmp_ss(a, b, 7) }
}
/// cmpps
#[inline]
pub fn _mm_cmpord_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.eq(a) & b.eq(b))
}
/// cmpss
#[inline]
pub fn _mm_cmpunord_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_cmp_ss(a, b, 3) }
}
/// cmpps
#[inline]
pub fn _mm_cmpunord_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.ne(a) | b.ne(b))
}
/// comiss
#[inline]
pub fn _mm_comieq_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_comieq_ss(a, b) }
}
/// comiss
#[inline]
pub fn _mm_comilt_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_comilt_ss(a, b) }
}
/// comiss
#[inline]
pub fn _mm_comile_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_comile_ss(a, b) }
}
/// comiss
#[inline]
pub fn _mm_comigt_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_comigt_ss(a, b) }
}
/// comiss
#[inline]
pub fn _mm_comige_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_comige_ss(a, b) }
}
/// comiss
#[inline]
pub fn _mm_comineq_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_comineq_ss(a, b) }
}
/// ucomiss
#[inline]
pub fn _mm_ucomieq_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_ucomieq_ss(a, b) }
}
/// ucomiss
#[inline]
pub fn _mm_ucomilt_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_ucomilt_ss(a, b) }
}
/// ucomiss
#[inline]
pub fn _mm_ucomile_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_ucomile_ss(a, b) }
}
/// ucomiss
#[inline]
pub fn _mm_ucomigt_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_ucomigt_ss(a, b) }
}
/// ucomiss
#[inline]
pub fn _mm_ucomige_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_ucomige_ss(a, b) }
}
/// ucomiss
#[inline]
pub fn _mm_ucomineq_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_ucomineq_ss(a, b) }
}
/// cvtss2si
#[inline]
pub fn _mm_cvtss_si32(a: __m128) -> i32 {
    unsafe { sse_cvtss2si(a) }
}
/// cvtss2si
#[inline]
pub fn _mm_cvt_ss2si(a: __m128) -> i32 {
    _mm_cvtss_si32(a)
}
/// cvtss2si
#[inline]
#[cfg(target_pointer_width = "64")]
pub fn _mm_cvtss_si64(a: __m128) -> i64 {
    unsafe { sse_cvtss2si64(a) }
}
/// cvttss2si
#[inline]
pub fn _mm_cvttss_si32(a: __m128) -> i32 {
    unsafe { sse_cvttss2si(a) }
}
/// cvttss2si
#[inline]
pub fn _mm_cvtt_ss2si(a: __m128) -> i32 {
    _mm_cvttss_si32(a)
}
/// cvttss2si
#[inline]
#[cfg(target_pointer_width = "64")]
pub fn _mm_cvttss_si64(a: __m128) -> i64 {
    unsafe { sse_cvttss2si64(a) }
}
/// cvtsi2ss
#[inline]
pub fn _mm_cvtsi32_ss(a: __m128, b: i32) -> __m128 {
    unsafe { sse_cvtsi2ss(a, b) }
}
/// cvtsi2ss
#[inline]
pub fn _mm_cvt_si2ss(a: __m128, b: i32) -> __m128 {
    _mm_cvtsi32_ss(a, b)
}
/// cvtsi2ss
#[inline]
#[cfg(target_pointer_width = "64")]
pub fn _mm_cvtsi64_ss(a: __m128, b: i64) -> __m128 {
    unsafe { sse_cvtsi642ss(a, b) }
}
#[inline]
pub fn _mm_cvtss_f32(a: __m128) -> f32 {
    a.extract(0)
}
/// movhps
#[inline]
pub unsafe fn _mm_loadh_pi(mut a: __m128, p: *const __m64) -> __m128 {
    // Cast to u8 to avoid assumptions about alignment.
    let p = p as *const u8;
    let pa = &mut a as *mut __m128 as *mut u8;
    copy_nonoverlapping(p, pa.offset(8), 8);
    a
}
/// movlps
#[inline]
pub unsafe fn _mm_loadl_pi(mut a: __m128, p: *const __m64) -> __m128 {
    // Cast to u8 to avoid assumptions about alignment.
    // FIXME: For some reason, this comes out as two instructions;
    // write this differently?
    let p = p as *const u8;
    let pa = &mut a as *mut __m128 as *mut u8;
    copy_nonoverlapping(p, pa, 8);
    a
}
/// movss
#[inline]
pub unsafe fn _mm_load_ss(p: *const f32) -> __m128 {
    let mut f: f32 = 0.;
    let p = p as *const u8;
    let pf = &mut f as *mut f32 as *mut u8;
    copy_nonoverlapping(p, pf, 4);
    __m128::new(f, 0., 0., 0.)
}
#[inline]
pub unsafe fn _mm_load1_ps(p: *const f32) -> __m128 {
    let mut f: f32 = 0.;
    let p = p as *const u8;
    let pf = &mut f as *mut f32 as *mut u8;
    copy_nonoverlapping(p, pf, 4);
    __m128::new(f, f, f, f)
}
#[inline]
pub unsafe fn _mm_load_ps1(p: *const f32) -> __m128 {
    _mm_load1_ps(p)
}
/// movaps
#[inline]
pub unsafe fn _mm_load_ps(p: *const f32) -> __m128 {
    let p = p as *mut __m128;
    *p
}
/// movups
#[inline]
pub unsafe fn _mm_loadu_ps(p: *const f32) -> __m128 {
    let mut a: __m128 = __m128::new(0., 0., 0., 0.);
    let p = p as *const u8;
    let pa = &mut a as *mut __m128 as *mut u8;
    copy_nonoverlapping(p, pa, 16);
    a
}
#[inline]
pub unsafe fn _mm_loadr_ps(p: *const f32) -> __m128 {
    let a = _mm_load_ps(p);
    __m128::new(a.extract(3), a.extract(2), a.extract(1), a.extract(0))
}
#[inline]
pub unsafe fn _mm_undefined_ps() -> __m128 {
    uninitialized()
}
#[inline]
pub fn _mm_set_ss(w: f32) -> __m128 {
    __m128::new(w, 0., 0., 0.)
}
#[inline]
pub fn _mm_set1_ps(w: f32) -> __m128 {
    __m128::new(w, w, w, w)
}
#[inline]
pub fn _mm_set_ps1(w: f32) -> __m128 {
    _mm_set1_ps(w)
}
#[inline]
pub fn _mm_set_ps(z: f32, y: f32, x: f32, w: f32) -> __m128 {
    __m128::new(w, x, y, z)
}
#[inline]
pub fn _mm_setr_ps(z: f32, y: f32, x: f32, w: f32) -> __m128 {
    __m128::new(z, y, x, w)
}
#[inline]
pub fn _mm_setzero_ps() -> __m128 {
    __m128::new(0., 0., 0., 0.)
}
/// movhps
#[inline]
pub unsafe fn _mm_storeh_pi(p: *mut __m64, a: __m128) {
    // Cast to u8 to avoid assumptions about alignment.
    // FIXME: For some reason, this comes out as two instructions;
    // write this differently?
    let p = p as *mut u8;
    let pa = &a as *const __m128 as *const u8;
    copy_nonoverlapping(pa.offset(8), p, 8)
}
/// movlps
#[inline]
pub unsafe fn _mm_storel_pi(p: *mut __m64, a: __m128) {
    // Cast to u8 to avoid assumptions about alignment.
    let p = p as *mut u8;
    let pa = &a as *const __m128 as *const u8;
    copy_nonoverlapping(pa, p, 8)
}
/// movss
#[inline]
pub unsafe fn _mm_store_ss(p: *mut f32, a: __m128) {
    // Cast to u8 to avoid assumptions about alignment.
    let p = p as *mut u8;
    let pa = &a as *const __m128 as *const u8;
    copy_nonoverlapping(pa, p, 4)
}
/// movups
#[inline]
pub unsafe fn _mm_storeu_ps(p: *mut f32, a: __m128) {
    // Cast to u8 to avoid assumptions about alignment.
    let p = p as *mut u8;
    let pa = &a as *const __m128 as *const u8;
    copy_nonoverlapping(pa, p, 16)
}
/// movaps
#[inline]
pub unsafe fn _mm_store_ps(p: *mut f32, a: __m128) {
    let p = p as *mut __m128;
    *p = a;
}
#[inline]
pub unsafe fn _mm_store1_ps(p: *mut f32, a: __m128) {
    let elem = a.extract(0);
    _mm_store_ps(p, __m128::new(elem, elem, elem, elem));
}
#[inline]
pub unsafe fn _mm_store_ps1(p: *mut f32, a: __m128) {
    _mm_store1_ps(p, a);
}
#[inline]
pub unsafe fn _mm_storer_ps(p: *mut f32, a: __m128) {
    _mm_store_ps(p, __m128::new(a.extract(3), a.extract(2), a.extract(1), a.extract(0)));
}
/// sfence
#[inline]
pub fn _mm_sfence() {
    unsafe { sse_sfence() }
}
/// stmxcsr
#[inline]
pub fn _mm_getcsr() -> u32 {
    unsafe {
        let mut result : u32 = uninitialized();
        sse_stmxcsr(&mut result as *mut u32 as *mut i8);
        result
    }
}
/// ldmxcsr
#[inline]
pub unsafe fn _mm_setcsr(i: u32) {
    sse_ldmxcsr(&i as *const u32 as *mut i8);
}
/// unpckhps
#[inline]
pub fn _mm_unpackhi_ps(a: __m128, b: __m128) -> __m128 {
    __m128::new(a.extract(2), b.extract(2), a.extract(3), b.extract(3))
}
/// unpcklps
#[inline]
pub fn _mm_unpacklo_ps(a: __m128, b: __m128) -> __m128 {
    __m128::new(a.extract(0), b.extract(0), a.extract(1), b.extract(1))
}
/// movss
#[inline]
pub fn _mm_move_ss(a: __m128, b: __m128) -> __m128 {
    a.replace(0, b.extract(0))
}
/// movhlps
#[inline]
pub fn _mm_movehl_ps(a: __m128, b: __m128) -> __m128 {
    // (This actually generates unpckhpd if SSE2 is available.)
    __m128::new(b.extract(2), b.extract(3), a.extract(2), a.extract(3))
}
/// movlhps
#[inline]
pub fn _mm_movelh_ps(a: __m128, b: __m128) -> __m128 {
    // (This actually generates unpcklpd if SSE2 is available.)
    __m128::new(a.extract(0), a.extract(1), b.extract(0), b.extract(1))
}
/// movmskps
#[inline]
pub fn _mm_movemask_ps(a: __m128) -> i32 {
    unsafe { sse_movmsk_ps(a) }
}
/// shufps
#[inline]
pub fn _mm_shuffle_ps(a: __m128, b: __m128, imm8: u32) -> __m128 {
    // FIXME: This is extremely ugly, but I don't know what the alternative is.
    let a1 = (imm8 >> 0) & 3;
    let a2 = (imm8 >> 2) & 3;
    let b1 = (imm8 >> 4) & 3;
    let b2 = (imm8 >> 6) & 3;
    __m128::new(a.extract(a1), a.extract(a2), b.extract(b1), b.extract(b2))
}
/// prefetchnta, prefetcht0, prefetcht1, prefetcht2
#[inline]
pub fn _mm_prefetch(p: *const i8, i: i32) {
    unsafe {
        match i {
            _MM_HINT_T0 => prefetch(p, 0, 3, 1),
            _MM_HINT_T1 => prefetch(p, 0, 2, 1),
            _MM_HINT_T2 => prefetch(p, 0, 1, 1),
            _ => prefetch(p, 0, 0, 1),
        }
    }
}
#[inline]
pub fn _MM_SHUFFLE(z: u32, y: u32, x: u32, w: u32) -> u32 {
    (z << 6) | (y << 4) | (x << 2) | w
}

/*
 TODO
unsigned int _MM_GET_EXCEPTION_MASK ()
unsigned int _MM_GET_EXCEPTION_STATE ()
unsigned int _MM_GET_FLUSH_ZERO_MODE ()
unsigned int _MM_GET_ROUNDING_MODE ()
void _MM_SET_EXCEPTION_MASK (unsigned int a)
void _MM_SET_EXCEPTION_STATE (unsigned int a)
void _MM_SET_FLUSH_ZERO_MODE (unsigned int a)
void _MM_SET_ROUNDING_MODE (unsigned int a)
_MM_TRANSPOSE4_PS (__m128 row0, __m128 row1, __m128 row2, __m128 row3)
*/

pub const _MM_EXCEPT_INVALID: u32 = 0x0001;
pub const _MM_EXCEPT_DENORM: u32 = 0x0002;
pub const _MM_EXCEPT_DIV_ZERO: u32 = 0x0004;
pub const _MM_EXCEPT_OVERFLOW: u32 = 0x0008;
pub const _MM_EXCEPT_UNDERFLOW: u32 = 0x0010;
pub const _MM_EXCEPT_INEXACT: u32 = 0x0020;
pub const _MM_EXCEPT_MASK: u32 = 0x003f;
pub const _MM_MASK_INVALID: u32 = 0x0080;
pub const _MM_MASK_DENORM: u32 = 0x0100;
pub const _MM_MASK_DIV_ZERO: u32 = 0x0200;
pub const _MM_MASK_OVERFLOW: u32 = 0x0400;
pub const _MM_MASK_UNDERFLOW: u32 = 0x0800;
pub const _MM_MASK_INEXACT: u32 = 0x1000;
pub const _MM_MASK_MASK: u32 = 0x1f80;
pub const _MM_ROUND_NEAREST: u32 = 0x0000;
pub const _MM_ROUND_DOWN: u32 = 0x2000;
pub const _MM_ROUND_UP: u32 = 0x4000;
pub const _MM_ROUND_TOWARD_ZERO: u32 = 0x6000;
pub const _MM_ROUND_MASK: u32 = 0x6000;
pub const _MM_FLUSH_ZERO_MASK: u32 = 0x8000;
pub const _MM_FLUSH_ZERO_ON: u32 = 0x8000;
pub const _MM_FLUSH_ZERO_OFF: u32 = 0x0000;

pub const _MM_HINT_T0 : i32 = 3;
pub const _MM_HINT_T1 : i32 = 2;
pub const _MM_HINT_T2 : i32 = 1;
pub const _MM_HINT_NTA : i32 = 0;

/// The methods in this module can't be implemented because Rust doesn't
/// expose nontemporal loads.
pub mod unimplemented_nontemporal {
    use __m128;
    /// movntps
    ///
    /// Not yet implemented.
    pub fn _mm_stream_ps(p: *mut f32, a: __m128) {
        unimplemented!()
    }
}

/// The methods in this module can't be implemented because Rust doesn't
/// expose the LLVM x86_mmx type.
pub mod unimplemented_mmx {
    use __m64;
    use __m128;

    /// cvtps2pi
    ///
    /// Not yet implemented.
    pub unsafe fn _mm_cvtps_pi32(a: __m128) -> __m64 {
        unimplemented!()
    }
    /// cvtps2pi
    ///
    /// Not yet implemented.
    pub unsafe fn _mm_cvt_ps2pi(a: __m128) -> __m64 {
        _mm_cvtps_pi32(a)
    }
    /// cvttps2pi
    ///
    /// Not yet implemented.
    pub unsafe fn _mm_cvttps_pi32(a: __m128) -> __m64 {
        unimplemented!()
    }
    /// cvttps2pi
    ///
    /// Not yet implemented.
    pub unsafe fn _mm_cvtt_ps2pi(a: __m128) -> __m64 {
        _mm_cvttps_pi32(a)
    }
    /// cvtpi2ps
    ///
    /// Not yet implemented.
    pub fn _mm_cvtpi32_ps(a: __m128, b: __m64) -> __m128 {
        unimplemented!()
    }
    /// cvtpi2ps
    ///
    /// Not yet implemented.
    pub fn _mm_cvt_pi2ps(a: __m128, b: __m64) -> __m128 {
        _mm_cvtpi32_ps(a, b)
    }
    /// movntq
    ///
    /// Not yet implemented.
    pub fn _mm_stream_pi(p: *mut __m64, a: __m64) {
        unimplemented!()
    }
    /// pmaxsw
    ///
    /// Not yet implemented.
    pub unsafe fn _mm_max_pi16(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// pmaxub
    ///
    /// Not yet implemented.
    pub unsafe fn _mm_max_pu8(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// pminsw
    ///
    /// Not yet implemented.
    pub unsafe fn _mm_min_pi16(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// pminub
    ///
    /// Not yet implemented.
    pub unsafe fn _mm_min_pu8(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// pmovmskb
    ///
    /// Not yet implemented.
    pub fn _mm_movemask_pi8(a: __m64) -> i32 {
        unimplemented!()
    }
    /// pmulhuw
    ///
    /// Not yet implemented.
    pub unsafe fn _mm_mulhi_pu16(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// maskmovq
    ///
    /// Not yet implemented.
    pub fn _mm_maskmove_si64(d: __m64, n: __m64, p: *mut i8) {
        unimplemented!()
    }
    /// pavgb
    ///
    /// Not yet implemented.
    pub unsafe fn _mm_avg_pu8(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// pavgw
    ///
    /// Not yet implemented.
    pub unsafe fn _mm_avg_pu16(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// psadbw
    ///
    /// Not yet implemented.
    pub unsafe fn _mm_sad_pu8(a: __m64, b: __m64) -> __m64 {
        unimplemented!()
    }
    /// Not yet implemented.
    pub fn _mm_cvtpi16_ps(a: __m64) -> __m128 {
        unimplemented!()
    }
    /// Not yet implemented.
    pub fn _mm_cvtpu16_ps(a: __m64) -> __m128 {
        unimplemented!()
    }
    /// Not yet implemented.
    pub fn _mm_cvtpi8_ps(a: __m64) -> __m128 {
        unimplemented!()
    }
    /// Not yet implemented.
    pub fn _mm_cvtpu8_ps(a: __m64) -> __m128 {
        unimplemented!()
    }
    /// Not yet implemented.
    pub fn _mm_cvtpi32x2_ps(a: __m64, b: __m64) -> __m128 {
        unimplemented!()
    }
    /// Not yet implemented.
    pub unsafe fn _mm_cvtps_pi16(a: __m128) -> __m64 {
        unimplemented!()
    }
    /// Not yet implemented.
    pub unsafe fn _mm_cvtps_pi8(a: __m128) -> __m64 {
        unimplemented!()
    }
    /// pextrw
    ///
    /// Not yet implemented.
    pub unsafe fn _mm_extract_pi16(a: __m64, imm8: i32) -> i32 {
        unimplemented!()
    }
    /// pinsrw
    ///
    /// Not yet implemented.
    pub unsafe fn _mm_insert_pi16(a: __m64, i: i32, imm8: i32) -> __m64 {
        unimplemented!()
    }
    /// pshufw
    ///
    /// Not yet implemented.
    pub unsafe fn _mm_shuffle_pi16(a: __m64, imm8: i32) -> __m64 {
        unimplemented!()
    }
    /// maskmovq
    ///
    /// Not yet implemented.
    pub unsafe fn _m_maskmovq(a: __m64, mask: __m64, mem_addr: *mut i8) {
        _mm_maskmove_si64(a, mask, mem_addr);
    }
    /// pavgb
    ///
    /// Not yet implemented.
    pub unsafe fn _m_pavgb(a: __m64, b: __m64) -> __m64 {
        _mm_avg_pu8(a, b)
    }
    /// pavgw
    ///
    /// Not yet implemented.
    pub unsafe fn _m_pavgw(a: __m64, b: __m64) -> __m64 {
        _mm_avg_pu16(a, b)
    }
    /// pextrw
    ///
    /// Not yet implemented.
    pub unsafe fn _m_pextrw(a: __m64, imm8: i32) -> i32 {
        _mm_extract_pi16(a, imm8)
    }
    /// pinsrw
    ///
    /// Not yet implemented.
    pub unsafe fn _m_pinsrw(a: __m64, i: i32, imm8: i32) -> __m64 {
        _mm_insert_pi16(a, i, imm8)
    }
    /// pmaxsw
    ///
    /// Not yet implemented.
    pub unsafe fn _m_pmaxsw(a: __m64, b: __m64) -> __m64 {
        _mm_max_pi16(a, b)
    }
    /// pmaxub
    ///
    /// Not yet implemented.
    pub unsafe fn _m_pmaxub(a: __m64, b: __m64) -> __m64 {
        _mm_max_pu8(a, b)
    }
    /// pminsw
    ///
    /// Not yet implemented.
    pub unsafe fn _m_pminsw(a: __m64, b: __m64) -> __m64 {
        _mm_min_pi16(a, b)
    }
    /// pminub
    ///
    /// Not yet implemented.
    pub unsafe fn _m_pminub(a: __m64, b: __m64) -> __m64 {
        _mm_min_pu8(a, b)
    }
    /// pmovmskb
    ///
    /// Not yet implemented.
    pub unsafe fn _m_pmovmskb(a: __m64) -> i32 {
        _mm_movemask_pi8(a)
    }
    /// pmulhuw
    ///
    /// Not yet implemented.
    pub unsafe fn _m_pmulhuw(a: __m64, b: __m64) -> __m64 {
        _mm_mulhi_pu16(a, b)
    }
    /// psadbw
    ///
    /// Not yet implemented.
    pub unsafe fn _m_psadbw(a: __m64, b: __m64) -> __m64 {
        _mm_sad_pu8(a, b)
    }
    /// pshufw
    ///
    /// Not yet implemented.
    pub unsafe fn _m_pshufw(a: __m64, imm8: i32) -> __m64 {
        _mm_shuffle_pi16(a, imm8)
    }
}
