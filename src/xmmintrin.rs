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

#[inline]
pub fn _mm_add_ss(a: __m128, b: __m128) -> __m128 {
    a.replace(0, a.extract(0) + b.extract(0))
}
#[inline]
pub fn _mm_add_ps(a: __m128, b: __m128) -> __m128 {
    a + b
}
#[inline]
pub fn _mm_sub_ss(a: __m128, b: __m128) -> __m128 {
    a.replace(0, a.extract(0) - b.extract(0))
}
#[inline]
pub fn _mm_sub_ps(a: __m128, b: __m128) -> __m128 {
    a - b
}
#[inline]
pub fn _mm_mul_ss(a: __m128, b: __m128) -> __m128 {
    a.replace(0, a.extract(0) * b.extract(0))
}
#[inline]
pub fn _mm_mul_ps(a: __m128, b: __m128) -> __m128 {
    a * b
}
#[inline]
pub fn _mm_div_ss(a: __m128, b: __m128) -> __m128 {
    a.replace(0, a.extract(0) / b.extract(0))
}
#[inline]
pub fn _mm_div_ps(a: __m128, b: __m128) -> __m128 {
    a / b
}
#[inline]
pub fn _mm_sqrt_ss(a: __m128) -> __m128 {
    unsafe { sse_sqrt_ss(a) }
}
#[inline]
pub fn _mm_sqrt_ps(a: __m128) -> __m128 {
    a.sqrt()
}
#[inline]
pub fn _mm_rcp_ss(a: __m128) -> __m128 {
    unsafe { sse_rcp_ss(a) }
}
#[inline]
pub fn _mm_rcp_ps(a: __m128) -> __m128 {
    a.approx_reciprocal()
}
#[inline]
pub fn _mm_rsqrt_ss(a: __m128) -> __m128 {
    unsafe { sse_rsqrt_ss(a) }
}
#[inline]
pub fn _mm_rsqrt_ps(a: __m128) -> __m128 {
    a.approx_rsqrt()
}
#[inline]
pub fn _mm_min_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_min_ss(a, b) }
}
#[inline]
pub fn _mm_min_ps(a: __m128, b: __m128) -> __m128 {
    a.min(b)
}
#[inline]
pub fn _mm_max_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_max_ss(a, b) }
}
#[inline]
pub fn _mm_max_ps(a: __m128, b: __m128) -> __m128 {
    a.max(b)
}
#[inline]
pub fn _mm_and_ps(a: __m128, b: __m128) -> __m128 {
    convert_m128i_to_m128(convert_m128_to_m128i(a) & convert_m128_to_m128i(b))
}
#[inline]
pub fn _mm_andnot_ps(a: __m128, b: __m128) -> __m128 {
    // FIXME: Not operator from simd doesn't get inlined?!
    convert_m128i_to_m128((convert_m128_to_m128i(a) ^ __m128i::splat(!0)) & convert_m128_to_m128i(b))
}
#[inline]
pub fn _mm_or_ps(a: __m128, b: __m128) -> __m128 {
    convert_m128i_to_m128(convert_m128_to_m128i(a) | convert_m128_to_m128i(b))
}
#[inline]
pub fn _mm_xor_ps(a: __m128, b: __m128) -> __m128 {
    convert_m128i_to_m128(convert_m128_to_m128i(a) ^ convert_m128_to_m128i(b))
}
#[inline]
pub fn _mm_cmpeq_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_cmp_ss(a, b, 0) }
}
#[inline]
pub fn _mm_cmpeq_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.eq(b))
}
#[inline]
pub fn _mm_cmplt_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_cmp_ss(a, b, 1) }
}
#[inline]
pub fn _mm_cmplt_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.lt(b))
}
#[inline]
pub fn _mm_cmple_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_cmp_ss(a, b, 2) }
}
#[inline]
pub fn _mm_cmple_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.le(b))
}
#[inline]
pub fn _mm_cmpgt_ss(a: __m128, b: __m128) -> __m128 {
    _mm_move_ss(a, _mm_cmplt_ss(b, a))
}
#[inline]
pub fn _mm_cmpgt_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.gt(b))
}
#[inline]
pub fn _mm_cmpge_ss(a: __m128, b: __m128) -> __m128 {
    _mm_move_ss(a, _mm_cmple_ss(b, a))
}
#[inline]
pub fn _mm_cmpge_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.ge(b))
}
#[inline]
pub fn _mm_cmpneq_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_cmp_ss(a, b, 4) }
}
#[inline]
pub fn _mm_cmpneq_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.ne(b))
}
#[inline]
pub fn _mm_cmpnlt_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_cmp_ss(a, b, 5) }
}
#[inline]
pub fn _mm_cmpnlt_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.ge(b) | (a.ne(a) | b.ne(b)))
}
#[inline]
pub fn _mm_cmpnle_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_cmp_ss(a, b, 6) }
}
#[inline]
pub fn _mm_cmpnle_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.gt(b) | (a.ne(a) | b.ne(b)))
}
#[inline]
pub fn _mm_cmpngt_ss(a: __m128, b: __m128) -> __m128 {
    _mm_move_ss(a, _mm_cmpnlt_ss(b, a))
}
#[inline]
pub fn _mm_cmpngt_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.le(b) | (a.ne(a) | b.ne(b)))
}
#[inline]
pub fn _mm_cmpnge_ss(a: __m128, b: __m128) -> __m128 {
    _mm_move_ss(a, _mm_cmpnle_ss(b, a))
}
#[inline]
pub fn _mm_cmpnge_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.lt(b) |(a.ne(a) | b.ne(b)))
}
#[inline]
pub fn _mm_cmpord_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_cmp_ss(a, b, 7) }
}
#[inline]
pub fn _mm_cmpord_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.eq(a) & b.eq(b))
}
#[inline]
pub fn _mm_cmpunord_ss(a: __m128, b: __m128) -> __m128 {
    unsafe { sse_cmp_ss(a, b, 3) }
}
#[inline]
pub fn _mm_cmpunord_ps(a: __m128, b: __m128) -> __m128 {
    convert_bool32fx4_to_m128(a.ne(a) | b.ne(b))
}
#[inline]
pub fn _mm_comieq_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_comieq_ss(a, b) }
}
#[inline]
pub fn _mm_comilt_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_comilt_ss(a, b) }
}
#[inline]
pub fn _mm_comile_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_comile_ss(a, b) }
}
#[inline]
pub fn _mm_comigt_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_comigt_ss(a, b) }
}
#[inline]
pub fn _mm_comige_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_comige_ss(a, b) }
}
#[inline]
pub fn _mm_comineq_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_comineq_ss(a, b) }
}
#[inline]
pub fn _mm_ucomieq_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_ucomieq_ss(a, b) }
}
#[inline]
pub fn _mm_ucomilt_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_ucomilt_ss(a, b) }
}
#[inline]
pub fn _mm_ucomile_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_ucomile_ss(a, b) }
}
#[inline]
pub fn _mm_ucomigt_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_ucomigt_ss(a, b) }
}
#[inline]
pub fn _mm_ucomige_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_ucomige_ss(a, b) }
}
#[inline]
pub fn _mm_ucomineq_ss(a: __m128, b: __m128) -> i32 {
    unsafe { sse_ucomineq_ss(a, b) }
}
#[inline]
pub fn _mm_cvtss_si32(a: __m128) -> i32 {
    unsafe { sse_cvtss2si(a) }
}
#[inline]
pub fn _mm_cvt_ss2si(a: __m128) -> i32 {
    _mm_cvtss_si32(a)
}
#[inline]
#[cfg(target_pointer_width = "64")]
pub fn _mm_cvtss_si64(a: __m128) -> i64 {
    unsafe { sse_cvtss2si64(a) }
}

#[inline]
pub fn _mm_cvttss_si32(a: __m128) -> i32 {
    unsafe { sse_cvttss2si(a) }
}
#[inline]
pub fn _mm_cvtt_ss2si(a: __m128) -> i32 {
    _mm_cvttss_si32(a)
}
#[inline]
#[cfg(target_pointer_width = "64")]
pub fn _mm_cvttss_si64(a: __m128) -> i64 {
    unsafe { sse_cvttss2si64(a) }
}
#[inline]
pub fn _mm_cvtsi32_ss(a: __m128, b: i32) -> __m128 {
    unsafe { sse_cvtsi2ss(a, b) }
}
#[inline]
pub fn _mm_cvt_si2ss(a: __m128, b: i32) -> __m128 {
    _mm_cvtsi32_ss(a, b)
}
#[inline]
pub fn _mm_cvtss_f32(a: __m128) -> f32 {
    a.extract(0)
}
#[inline]
pub unsafe fn _mm_loadh_pi(mut a: __m128, p: *const __m64) -> __m128 {
    // Cast to u8 to avoid assumptions about alignment.
    let p = p as *const u8;
    let pa = &mut a as *mut __m128 as *mut u8;
    copy_nonoverlapping(p, pa.offset(8), 8);
    a
}
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
#[inline]
pub unsafe fn _mm_load_ps(p: *const f32) -> __m128 {
    let p = p as *mut __m128;
    *p
}
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
#[inline]
pub unsafe fn _mm_storeh_pi(p: *mut __m64, a: __m128) {
    // Cast to u8 to avoid assumptions about alignment.
    // FIXME: For some reason, this comes out as two instructions;
    // write this differently?
    let p = p as *mut u8;
    let pa = &a as *const __m128 as *const u8;
    copy_nonoverlapping(pa.offset(8), p, 8)
}
#[inline]
pub unsafe fn _mm_storel_pi(p: *mut __m64, a: __m128) {
    // Cast to u8 to avoid assumptions about alignment.
    let p = p as *mut u8;
    let pa = &a as *const __m128 as *const u8;
    copy_nonoverlapping(pa, p, 8)
}
#[inline]
pub unsafe fn _mm_store_ss(p: *mut f32, a: __m128) {
    // Cast to u8 to avoid assumptions about alignment.
    let p = p as *mut u8;
    let pa = &a as *const __m128 as *const u8;
    copy_nonoverlapping(pa, p, 4)
}
#[inline]
pub unsafe fn _mm_storeu_ps(p: *mut f32, a: __m128) {
    // Cast to u8 to avoid assumptions about alignment.
    let p = p as *mut u8;
    let pa = &a as *const __m128 as *const u8;
    copy_nonoverlapping(pa, p, 16)
}
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
#[inline]
pub fn _mm_sfence() {
    unsafe { sse_sfence() }
}
#[inline]
pub fn _mm_getcsr() -> u32 {
    unsafe {
        let mut result : u32 = uninitialized();
        sse_stmxcsr(&mut result as *mut u32 as *mut i8);
        result
    }
}
#[inline]
pub unsafe fn _mm_setcsr(i: u32) {
    sse_ldmxcsr(&i as *const u32 as *mut i8);
}
#[inline]
pub fn _mm_unpackhi_ps(a: __m128, b: __m128) -> __m128 {
    __m128::new(a.extract(2), b.extract(2), a.extract(3), b.extract(3))
}
#[inline]
pub fn _mm_unpacklo_ps(a: __m128, b: __m128) -> __m128 {
    __m128::new(a.extract(0), b.extract(0), a.extract(1), b.extract(1))
}
#[inline]
pub fn _mm_move_ss(a: __m128, b: __m128) -> __m128 {
    a.replace(0, b.extract(0))
}
#[inline]
pub fn _mm_movehl_ps(a: __m128, b: __m128) -> __m128 {
    //(This actually generates unpckhpd if SSE2 is available.)
    __m128::new(b.extract(2), b.extract(3), a.extract(2), a.extract(3))
}
#[inline]
pub fn _mm_movelh_ps(a: __m128, b: __m128) -> __m128 {
    //(This actually generates unpcklpd if SSE2 is available.)
    __m128::new(a.extract(0), a.extract(1), b.extract(0), b.extract(1))
}
#[inline]
pub fn _mm_movemask_ps(a: __m128) -> i32 {
    unsafe { sse_movmsk_ps(a) }
}
#[inline]
pub fn _mm_shuffle_ps(a: __m128, b: __m128, imm8: u32) -> __m128 {
    // FIXME: This is extremely ugly, but I don't know what the alternative is.
    let a1 = (imm8 >> 0) & 3;
    let a2 = (imm8 >> 2) & 3;
    let b1 = (imm8 >> 4) & 3;
    let b2 = (imm8 >> 6) & 3;
    __m128::new(a.extract(a1), a.extract(a2), b.extract(b1), b.extract(b2))
}
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

#define         _MM_SHUFFLE(z, y, x, w)   (((z) << 6) | ((y) << 4) | ((x) << 2) | (w))
#define         _MM_EXCEPT_INVALID   (0x0001)
#define         _MM_EXCEPT_DENORM   (0x0002)
#define         _MM_EXCEPT_DIV_ZERO   (0x0004)
#define         _MM_EXCEPT_OVERFLOW   (0x0008)
#define         _MM_EXCEPT_UNDERFLOW   (0x0010)
#define         _MM_EXCEPT_INEXACT   (0x0020)
#define         _MM_EXCEPT_MASK   (0x003f)
#define         _MM_MASK_INVALID   (0x0080)
#define         _MM_MASK_DENORM   (0x0100)
#define         _MM_MASK_DIV_ZERO   (0x0200)
#define         _MM_MASK_OVERFLOW   (0x0400)
#define         _MM_MASK_UNDERFLOW   (0x0800)
#define         _MM_MASK_INEXACT   (0x1000)
#define         _MM_MASK_MASK   (0x1f80)
#define         _MM_ROUND_NEAREST   (0x0000)
#define         _MM_ROUND_DOWN   (0x2000)
#define         _MM_ROUND_UP   (0x4000)
#define         _MM_ROUND_TOWARD_ZERO   (0x6000)
#define         _MM_ROUND_MASK   (0x6000)
#define         _MM_FLUSH_ZERO_MASK   (0x8000)
#define         _MM_FLUSH_ZERO_ON   (0x8000)
#define         _MM_FLUSH_ZERO_OFF   (0x0000)
*/

pub const _MM_HINT_T0 : i32 = 3;
pub const _MM_HINT_T1 : i32 = 2;
pub const _MM_HINT_T2 : i32 = 1;
pub const _MM_HINT_NTA : i32 = 0;

/// The methods in this module can't be implemented because Rust doesn't
/// expose nontemporal loads.
pub mod unimplemented_nontemporal {
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
/// Not yet implemented.
pub unsafe fn _mm_cvtps_pi32(a: __m128) -> __m64 {
    unimplemented!()
}
/// Not yet implemented.
pub unsafe fn _mm_cvt_ps2pi(a: __m128) -> __m64 {
    unimplemented!()
}
/// Not yet implemented.
pub unsafe fn _mm_cvttps_pi32(a: __m128) -> __m64 {
    unimplemented!()
}
/// Not yet implemented.
pub unsafe fn _mm_cvtt_ps2pi(a: __m128) -> __m64 {
    unimplemented!()
}
/// Not yet implemented.
pub fn _mm_cvtpi32_ps(a: __m128, b: __m64) -> __m128 {
    unimplemented!()
}
/// Not yet implemented.
pub fn _mm_cvt_pi2ps(a: __m128, b: __m64) -> __m128 {
    unimplemented!()
}
/// Not yet implemented.
pub fn _mm_stream_pi(p: *mut __m64, a: __m64) {
    unimplemented!()
}
/// Not yet implemented.
pub unsafe fn _mm_max_pi16(a: __m64, b: __m64) -> __m64 {
    unimplemented!()
}
/// Not yet implemented.
pub unsafe fn _mm_max_pu8(a: __m64, b: __m64) -> __m64 {
    unimplemented!()
}
/// Not yet implemented.
pub unsafe fn _mm_min_pi16(a: __m64, b: __m64) -> __m64 {
    unimplemented!()
}
/// Not yet implemented.
pub unsafe fn _mm_min_pu8(a: __m64, b: __m64) -> __m64 {
    unimplemented!()
}
/// Not yet implemented.
pub fn _mm_movemask_pi8(a: __m64) -> i32 {
    unimplemented!()
}
/// Not yet implemented.
pub unsafe fn _mm_mulhi_pu16(a: __m64, b: __m64) -> __m64 {
    unimplemented!()
}
/// Not yet implemented.
pub fn _mm_maskmove_si64(d: __m64, n: __m64, p: *mut i8) {
    unimplemented!()
}
/// Not yet implemented.
pub unsafe fn _mm_avg_pu8(a: __m64, b: __m64) -> __m64 {
    unimplemented!()
}
/// Not yet implemented.
pub unsafe fn _mm_avg_pu16(a: __m64, b: __m64) -> __m64 {
    unimplemented!()
}
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
/// Not yet implemented.
pub unsafe fn _mm_extract_pi16(a: __m64, imm8: i32) -> i32 {
    unimplemented!()
}
/// Not yet implemented.
pub unsafe fn _mm_insert_pi16(a: __m64, i: i32, imm8: i32) -> __m64 {
    unimplemented!()
}
/// Not yet implemented.
pub unsafe fn _mm_shuffle_pi16(a: __m64, imm8: i32) -> __m64 {
    unimplemented!()
}
/// Not yet implemented.
pub unsafe fn _m_maskmovq(a: __m64, mask: __m64, mem_addr: *mut i8) {
    _mm_maskmove_si64(a, mask, mem_addr);
}
/// Not yet implemented.
pub unsafe fn _m_pavgb(a: __m64, b: __m64) -> __m64 {
    _mm_avg_pu8(a, b)
}
/// Not yet implemented.
pub unsafe fn _m_pavgw(a: __m64, b: __m64) -> __m64 {
    _mm_avg_pu16(a, b)
}
/// Not yet implemented.
pub unsafe fn _m_pextrw(a: __m64, imm8: i32) -> i32 {
    _mm_extract_pi16(a, imm8)
}
/// Not yet implemented.
pub unsafe fn _m_pinsrw (a: __m64, i: i32, imm8: i32) -> __m64 {
    _mm_insert_pi16(a, i, imm8)
}
/// Not yet implemented.
pub unsafe fn _m_pmaxsw(a: __m64, b: __m64) -> __m64 {
    _mm_max_pi16(a, b)
}
/// Not yet implemented.
pub unsafe fn _m_pmaxub(a: __m64, b: __m64) -> __m64 {
    _mm_max_pu8(a, b)
}
/// Not yet implemented.
pub unsafe fn _m_pminsw(a: __m64, b: __m64) -> __m64 {
    _mm_min_pi16(a, b)
}
/// Not yet implemented.
pub unsafe fn _m_pminub(a: __m64, b: __m64) -> __m64 {
    _mm_min_pu8(a, b)
}
/// Not yet implemented.
pub unsafe fn _m_pmovmskb(a: __m64) -> i32 {
    _mm_movemask_pi8(a)
}
/// Not yet implemented.
pub unsafe fn _m_pmulhuw(a: __m64, b: __m64) -> __m64 {
    _mm_mulhi_pu16(a, b)
}
/// Not yet implemented.
pub unsafe fn _m_psadbw(a: __m64, b: __m64) -> __m64 {
    _mm_sad_pu8(a, b)
}
/// Not yet implemented.
pub unsafe fn _m_pshufw(a: __m64, imm8: i32) -> __m64 {
    _mm_shuffle_pi16(a, imm8)
}
}
