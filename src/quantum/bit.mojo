# For now, Bits and Qubits are just Ints

@register_passable('trivial')
struct Clbit(StringableCollectionElement, KeyElement, Formattable):
    var idx: Int

    @always_inline
    fn __init__(inout self, idx: Int):
        self.idx = idx
    
    @no_inline
    fn __str__(self) -> String:
        return 'c' + str(self.idx)

    @no_inline
    fn format_to(self, inout writer: Formatter):
        writer.write(self.__str__())

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
struct Qubit(StringableCollectionElement, KeyElement, Formattable):
    var idx: Int

    @always_inline
    fn __init__(inout self, idx: Int):
        self.idx = idx
    
    @no_inline
    fn __str__(self) -> String:
        return 'q' + str(self.idx)

    @no_inline
    fn format_to(self, inout writer: Formatter):
        writer.write(self.__str__())

    @always_inline
    fn __hash__(self) -> UInt:
        return hash(self.idx)

    @always_inline
    fn __eq__(self, other: Self) -> Bool:
        return self.idx == other.idx

    @always_inline
    fn __ne__(self, other: Self) -> Bool:
        return self.idx != other.idx