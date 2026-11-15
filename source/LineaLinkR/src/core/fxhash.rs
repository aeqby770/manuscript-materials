use std::hash::{BuildHasherDefault, Hasher};

/// Minimal FxHasher — fast multiply-shift hash for integer keys.
/// Equivalent throughput to `rustc-hash` without the network dependency.
#[derive(Default, Clone)]
pub struct FxHasher {
    hash: u64,
}

const K: u64 = 0x517c_c1b7_2722_0a95;

impl FxHasher {
    #[inline]
    fn add(&mut self, v: u64) {
        // hash = hash.rotate_left(5) ^ i) * K
        self.hash = (self.hash.rotate_left(5) ^ v).wrapping_mul(K);
    }
}

impl Hasher for FxHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.hash
    }

    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        for chunk in bytes.chunks_exact(8) {
            self.add(u64::from_ne_bytes(chunk.try_into().unwrap()));
        }
        let rem = bytes.len() % 8;
        if rem > 0 {
            let tail = bytes[bytes.len() - rem..]
                .iter()
                .fold(0u64, |acc, &b| (acc << 8) | b as u64);
            self.add(tail);
        }
    }

    #[inline]
    fn write_u8(&mut self, v: u8) {
        self.add(v as u64);
    }
    #[inline]
    fn write_u16(&mut self, v: u16) {
        self.add(v as u64);
    }
    #[inline]
    fn write_u32(&mut self, v: u32) {
        self.add(v as u64);
    }
    #[inline]
    fn write_u64(&mut self, v: u64) {
        self.add(v);
    }
    #[inline]
    fn write_usize(&mut self, v: usize) {
        self.add(v as u64);
    }
    #[inline]
    fn write_i8(&mut self, v: i8) {
        self.add(v as u64);
    }
    #[inline]
    fn write_i16(&mut self, v: i16) {
        self.add(v as u64);
    }
    #[inline]
    fn write_i32(&mut self, v: i32) {
        self.add(v as u64);
    }
    #[inline]
    fn write_i64(&mut self, v: i64) {
        self.add(v as u64);
    }
    #[inline]
    fn write_isize(&mut self, v: isize) {
        self.add(v as u64);
    }
}

pub type FxHashMap<K, V> = std::collections::HashMap<K, V, BuildHasherDefault<FxHasher>>;
pub type FxHashSet<T> = std::collections::HashSet<T, BuildHasherDefault<FxHasher>>;
