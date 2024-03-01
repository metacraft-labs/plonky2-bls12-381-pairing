pub const BLS_X: u64 = 0xd201_0000_0001_0000;
pub const BLS_X_IS_NEGATIVE: bool = true;

// Canonical signed digit decomposition (Non-Adjacent form) of 6x₀+2 = 29793968203157093288  little endian
pub const NAF_DIGIT: [u8; 66] = [
    0, 0, 0, 1, 0, 1, 0, 2, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 2, 0,
    0, 1, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 1, 0, 2, 0, 0, 0, 2, 0, 2, 0, 0, 0, 1, 0, 2,
    0, 1,
];
