# immintrin

This crate implements the names defined by Intel for SIMD and related instructions on x86
(functions with an `_mm_` prefix). These are useful for porting code, and because people
are familiar with them.  See https://github.com/rust-lang/rfcs/issues/1639 for discussion.

Currently only SSE1 intrinsics are implemented; this crate is a work in progress.
Requires a nightly compiler for SIMD and various unstable intrinsics.

Licensed under the Apache License, Version 2.0
<LICENSE-APACHE or
http://www.apache.org/licenses/LICENSE-2.0> or the MIT
license <LICENSE-MIT or http://opensource.org/licenses/MIT>,
at your option. All files in the project carrying such
notice may not be copied, modified, or distributed except
according to those terms.
