# For now, Bits and Qubits are just Ints

@register_passable('trivial')
struct Clbit(Stringable, CollectionElement, KeyElement, Writable):
    var idx: Int

    @always_inline
    @implicit
    fn __init__(out self, idx: Int):
        self.idx = idx
    
    @no_inline
    fn __str__(self) -> String:
        return 'c' + String(self.idx)

    @no_inline
    fn write_to[W: Writer](self, mut writer: W):
        writer.write(String(self))

    @always_inline
    fn __hash__(self) -> UInt:
        return hash(self.idx)

    @always_inline
    fn __eq__(self, other: Self) -> Bool:
        return self.idx == other.idx

    @always_inline
    fn __ne__(self, other: Self) -> Bool:
        return self.idx != other.idx

@register_passable('trivial')
struct Qubit(Stringable, CollectionElement, KeyElement, Writable):
    var idx: Int

    @always_inline
    @implicit
    fn __init__(out self, idx: Int):
        self.idx = idx
    
    @no_inline
    fn __str__(self) -> String:
        return 'q' + String(self.idx)

    @no_inline
    fn write_to[W: Writer](self, mut writer: W):
        writer.write(String(self))

    @always_inline
    fn __hash__(self) -> UInt:
        return hash(self.idx)

    @always_inline
    fn __eq__(self, other: Self) -> Bool:
        return self.idx == other.idx

    @always_inline
    fn __ne__(self, other: Self) -> Bool:
        return self.idx != other.idx