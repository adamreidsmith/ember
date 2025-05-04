import math

from ..config import DEFAULT_TOL

alias ComplexScalar = ComplexSIMD[size=1]

alias ComplexSIMD16 = ComplexSIMD[DType.float16]
alias ComplexSIMD32 = ComplexSIMD[DType.float32]
alias ComplexSIMD64 = ComplexSIMD[DType.float64]
alias BComplexSIMD16 = ComplexSIMD[DType.bfloat16]

alias Complex16 = ComplexSIMD16[1]
alias Complex32 = ComplexSIMD32[1]
alias Complex64 = ComplexSIMD64[1]
alias BComplex16 = BComplexSIMD16[1]

alias _halfpi =  math.pi / 2


@register_passable('trivial')
struct ComplexSIMD[type: DType, size: Int](
        Defaultable, 
        Powable, 
        Roundable, 
        Absable,
        Intable,
        Stringable,
        RepresentableCollectionElement,
        Hashable,
        Sized,
        Writable,
    ):
    alias Coef = SIMD[Self.type, Self.size]
    alias Lane = ComplexSIMD[Self.type, 1]

    var re: Self.Coef
    var im: Self.Coef

    # Initialize ######################

    @always_inline
    fn __init__(out self):
        '''Initialize a complex number to 0 + 0i.'''
        self.re = 0
        self.im = 0        
    
    @always_inline
    @implicit
    fn __init__(out self, re: Self.Coef):
        '''Initialize a complex number with SIMD coefficients.'''
        self.re = re
        self.im = 0
    
    @always_inline
    fn __init__(out self, re: Self.Coef, im: Self.Coef):
        '''Initialize a complex number with SIMD coefficients.'''
        self.re = re
        self.im = im

    @always_inline
    fn __init__[__: None = None](out self, re: Scalar[Self.type], im: Scalar[Self.type]):
        '''Initialize a complex number with SIMD coefficients.'''
        self.re = re
        self.im = im
    
    # @always_inline
    # fn __init__(out self, t: (FloatLiteral, FloatLiteral)):
    #     '''Initialize a complex number with float literals.'''
    #     self.re = t[0]
    #     self.im = t[1]
    
    @always_inline
    fn __init__(out self, t: (Float16, Float16)):
        '''Initialize a complex number with float literals.'''
        self.re = t[0].cast[Self.type]()
        self.im = t[1].cast[Self.type]()
            
    @always_inline
    fn __init__(out self, t: (Float32, Float32)):
        '''Initialize a complex number with float literals.'''
        self.re = t[0].cast[Self.type]()
        self.im = t[1].cast[Self.type]()
            
    @always_inline
    @implicit
    fn __init__(out self, t: (Float64, Float64)):
        '''Initialize a complex number with float literals.'''
        self.re = t[0].cast[Self.type]()
        self.im = t[1].cast[Self.type]()

    @always_inline
    @implicit
    fn __init__(out self, t: (Int, Int)):
        '''Initialize a complex number with integer values.'''
        self.re = t[0]
        self.im = t[1]
    
    # @always_inline
    # fn __init__(out self, t: (IntLiteral, IntLiteral)):
    #     '''Initialize a complex number with integer literals.'''
    #     self.re = t[0]
    #     self.im = t[1]
    
    @always_inline
    fn __init__[__: None = None](out self, re: FloatLiteral):
        '''Initialize a complex number with a real float literal.'''
        self.re = re
        self.im = 0
    
    @always_inline
    fn __init__[__: None = None](out self, re: Float16):
        '''Initialize a complex number with a real float literal.'''
        self.re = re.cast[Self.type]()
        self.im = 0

    @always_inline
    fn __init__[__: None = None](out self, re: Float32):
        '''Initialize a complex number with a real float literal.'''
        self.re = re.cast[Self.type]()
        self.im = 0

    @always_inline
    fn __init__[__: None = None](out self, re: Float64):
        '''Initialize a complex number with a real float literal.'''
        self.re = re.cast[Self.type]()
        self.im = 0

    @always_inline
    @implicit
    fn __init__[__: None = None](out self, re: Int):
        '''Initialize a complex number with a real integer value.'''
        self.re = re
        self.im = 0
    
    # @always_inline
    # fn __init__[__: None = None](out self, re: IntLiteral):
    #     '''Initialize a complex number with a real integer literal.'''
    #     self.re = re
    #     self.im = 0

    # Static construction #############

    @always_inline
    @staticmethod
    fn from_polar(r: Self.Coef, theta: Self.Coef) -> Self:
        '''Create a complex number from polar coordinates (r, θ).'''
        return Self(r * math.cos(theta), r * math.sin(theta))
    
    @always_inline
    @staticmethod
    fn from_exp(r: Self.Coef, theta: Self.Coef) -> Self:
        '''Create a complex number from exponential form r*e^(iθ).'''
        return Self.from_polar(r, theta)
    
    @always_inline
    @staticmethod
    fn zero() -> Self:
        '''Create a complex number representing zero (0 + 0i).'''
        return Self(0, 0)
    
    @always_inline
    @staticmethod
    fn one() -> Self:
        '''Create a complex number representing one (1 + 0i).'''
        return Self(1, 0)
    
    @always_inline
    @staticmethod
    fn i() -> Self:
        '''Create the imaginary unit i (0 + 1i).'''
        return Self(0, 1)
    
    # Format ##########################

    @no_inline
    @staticmethod
    fn _fmt_float(owned st: String, max_digits: Int) -> String:
        if 'nan' in st or 'inf' in st: 
            return st
        var sign: String = ''
        if not st[0].isdigit():
            sign = '-' if st[0] == '-' else ''
            s = st[1:]
        else:
            s = st
        if 'e' in s:
            var e_idx: Int = s.find('e') 
            var mantissa: String = s[2:e_idx]
            while mantissa and (len(mantissa) > max_digits - 1 or mantissa[-1] == '0'):
                mantissa = mantissa[:-1]
            if not mantissa:
                mantissa = '0'
            return sign + s[:2] + mantissa[:max_digits - 1] + s[e_idx:]
        if '.' in s:
            var dot_idx: Int = s.find('.')
            var int_part: String = s[:dot_idx]
            var dec_part: String = s[dot_idx + 1:]
            if len(int_part) > max_digits:
                return (
                    sign + int_part[0] + '.' + int_part[1:max_digits] 
                    + 'e+' + String(len(int_part) - 1).rjust(2, '0')
                )
            if len(int_part) == max_digits:
                return sign + int_part
            return sign + int_part + '.' + dec_part[:max_digits - len(int_part)]
        if len(s) > max_digits:
            return sign + s[0] + '.' + s[1:max_digits] + 'e+' + String(len(s) - 1).rjust(2, '0')
        return sign + s

    @no_inline
    fn _str1(self, z: ComplexScalar[Self.type], max_digits: Int) -> String:
        if max_digits > 0:
            var re_str: String = self._fmt_float(String(z.re), max_digits)
            var im_str: String = self._fmt_float(String(z.im), max_digits) + 'i'
            if not im_str.startswith('-'):
                re_str += '+'
            return re_str + im_str
        var re_str: String = String(z.re)
        var im_str: String = String(z.im) + 'i'
        if not im_str.startswith('-'):
            re_str += '+'
        return re_str + im_str 

    @always_inline
    fn __str__(self) -> String:
        '''Formats the Complex as a String.'''
        return self.__str__(max_digits=8)

    @no_inline
    fn __str__[__: None = None](self, max_digits: Int) -> String:
        if Self.size == 1:
            return self._str1(self[0], max_digits)
        
        var out: String = '['
        for idx in range(Self.size):
            if idx > 0:
                out += ' '
            out += self._str1(self[idx], max_digits) + '\n'
        out = out[:-1]
        out += ']'
        return out
    
    @always_inline
    fn __repr__(self) -> String:
        '''Formats the Complex as a String.'''
        return self.__str__(max_digits=0)
    
    @no_inline
    fn write_to[W: Writer](self, mut writer: W):
        '''Formats the Complex as a String.'''
        writer.write(String(self))

    # Arithmetic ######################

    @always_inline
    fn __neg__(self) -> Self:
        '''Defines the `-` negative unary operator. Returns -self.'''
        return Self(-self.re, -self.im)

    @always_inline
    fn __pos__(self) -> Self:
        '''Defines the `+` positive unary operator. Returns +self.'''
        return Self(self.re, self.im)
    
    @always_inline
    fn __add__(self, other: Self) -> Self:
        '''Defines the `+` plus operator. Returns self + other.'''
        return Self(self.re + other.re, self.im + other.im)
    
    @always_inline
    fn __add__[__: None = None](self, other: Self.Lane) -> Self:
        '''Defines the `+` plus operator. Returns self + other.'''
        return Self(self.re + other.re, self.im + other.im)
    
    @always_inline
    fn __add__(self, other: Self.Coef) -> Self:
        '''Defines the `+` plus operator. Returns self + other.'''
        return Self(self.re + other, self.im)
    
    @always_inline
    fn __sub__(self, other: Self) -> Self:
        '''Defines the `-` minus operator. Returns self - other.'''
        return Self(self.re - other.re, self.im - other.im)

    @always_inline
    fn __sub__[__: None = None](self, other: Self.Lane) -> Self:
        '''Defines the `-` minus operator. Returns self - other.'''
        return Self(self.re - other.re, self.im - other.im)

    @always_inline
    fn __sub__[__: None = None](self, other: Self.Coef) -> Self:
        '''Defines the `-` minus operator. Returns self - other.'''
        return Self(self.re - other, self.im)
    
    @always_inline
    fn __mul__(self, other: Self) -> Self:
        '''Defines the `*` product operator. Returns self * other.'''
        return Self(
            self.re * other.re - self.im * other.im, self.re * other.im + self.im * other.re
        )

    @always_inline
    fn __mul__[__: None = None](self, other: Self.Lane) -> Self:
        '''Defines the `*` product operator. Returns self * other.'''
        return Self(
            self.re * other.re - self.im * other.im, self.re * other.im + self.im * other.re
        )

    @always_inline
    fn __mul__(self, other: Self.Coef) -> Self:
        '''Defines the `*` product operator. Returns self * other.'''
        return Self(self.re * other, self.im * other)
    
    @always_inline
    fn __truediv__(self, other: Self) -> Self:
        '''Defines the `/` divide operator. Returns self / other.'''
        var denom: Self.Coef = 1 / other.squared_norm()
        return Self(
            (self.re * other.re + self.im * other.im) * denom, 
            (self.im * other.re - self.re * other.im) * denom,
        )

    @always_inline
    fn __truediv__[__: None = None](self, other: Self.Lane) -> Self:
        '''Defines the `/` divide operator. Returns self / other.'''
        var denom: Self.Coef = 1 / other.squared_norm()
        return Self(
            (self.re * other.re + self.im * other.im) * denom, 
            (self.im * other.re - self.re * other.im) * denom,
        )
    
    @always_inline
    fn __truediv__(self, other: Self.Coef) -> Self:
        '''Defines the `/` divide operator. Returns self / other.'''
        return Self(self.re / other, self.im / other)
    
    @always_inline
    fn __floordiv__(self, other: Self) -> Self:
        '''Defines the `//` floor divide operator. Returns round(q.real()) + round(q.imag()) * i 
        for q = self / other.
        '''
        return round(self / other)

    @always_inline
    fn __floordiv__[__: None = None](self, other: Self.Lane) -> Self:
        '''Defines the `//` floor divide operator. Returns round(q.real()) + round(q.imag()) * i
        for q = self / other.
        '''
        return round(self / other)
    
    @always_inline
    fn __floordiv__(self, other: Self.Coef) -> Self:
        '''Defines the `//` floor divide operator. Returns round(q.real()) + round(q.imag()) * i
        for q = self / other.
        '''
        return round(self / other)

    @always_inline
    fn __mod__(self, other: Self) -> Self:
        '''Defines the `%` modulo operator. Returns self - (self // other) * other.'''
        return self - (self // other) * other
    
    @always_inline
    fn __mod__[__: None = None](self, other: Self.Lane) -> Self:
        '''Defines the `%` modulo operator. Returns self - (self // other) * other.'''
        return self - (self // other) * other
    
    @always_inline
    fn __mod__(self, other: Self.Coef) -> Self:
        '''Defines the `%` modulo operator. Returns self - (self // other) * other.'''
        return self - (self // other) * other
    
    @always_inline
    fn __divmod__(self, other: Self) -> Tuple[Self, Self]:
        '''Defines the divmod operator. Returns the tuple (self // other, self % other).'''
        return (self // other, self % other)
    
    @always_inline
    fn __divmod__[__: None = None](self, other: Self.Lane) -> Tuple[Self, Self]:
        '''Defines the divmod operator. Returns the tuple (self // other, self % other).'''
        return (self // other, self % other)
    
    @always_inline
    fn __divmod__(self, other: Self.Coef) -> Tuple[Self, Self]:
        '''Defines the divmod operator. Returns the tuple (self // other, self % other).'''
        return (self // other, self % other)

    @always_inline
    fn __abs__(self) -> Self:
        '''Returnes the absolute value of the complex number as a complex number
        with no imaginary part.
        '''
        return Self(math.sqrt(self.re * self.re + self.im * self.im))
    
    @always_inline
    fn __round__(self) -> Self:
        '''Round the real and imaginary parts of the complex number.'''
        return Self(round(self.re), round(self.im))

    @always_inline
    fn __round__(self, ndigits: Int) -> Self:
        '''Round the real and imaginary parts of the complex number.'''
        return Self(round(self.re, ndigits), round(self.im, ndigits))

    fn __pow__(self, other: Self) -> Self:
        '''Defines the `**` power operator. Returns self ** other.'''
        if other.im == 0:
            if other.re == 0:
                return Self(1, 0)
            elif other.re == 1:
                return self
            elif other.re == 2:
                return self * self
            elif other.re == 3:
                return self * self * self
        if self.re == 0 and self.im == 0:
            return Self(0, 0)
        return (other * self.log()).exp()
    
    @always_inline
    fn __pow__[__: None = None](self, other: Self.Lane) -> Self:
        '''Defines the `**` power operator. Returns self ** other.'''
        return self ** Self(other.re, other.im)
    
    @always_inline
    fn __pow__(self, other: Self.Coef) -> Self:
        '''Defines the `**` power operator. Returns self ** other.'''
        return self ** Self(other)
    
    @staticmethod
    fn _fast_int_pow(owned a: Self, b: Int) -> Self:
        if b == 1:
            return a
        elif b == 2:
            return a * a
        var result: Self = Self(1, 0)
        var p: Int = b
        while p:
            if p % 2 == 1:
                result *= a
            p //= 2
            a *= a
        return result
    
    @always_inline
    fn __pow__(self, b: Int) -> Self:
        '''Defines the `**` power operator. Returns self ** b.'''
        if b == 0:
            return Self(1, 0)
        if b < 0:
            return Self._fast_int_pow(self.copy(), -b).reciprocal()
        return Self._fast_int_pow(self.copy(), b)
    
    @always_inline
    fn __pow__(self, b: IntLiteral) -> Self:
        '''Defines the `**` power operator. Returns self ** b.'''
        return self ** Int(b)
    
    # Reverse arithmetic ##############

    @always_inline
    fn __radd__(self, other: Self) -> Self:
        '''Defines the reverse `+` plus operator. Returns other + self.'''
        return other + self

    @always_inline
    fn __radd__[__: None = None](self, other: Self.Lane) -> Self:
        '''Defines the reverse `+` plus operator. Returns other + self.'''
        return Self(other.re + self.re, other.im + self.im)

    @always_inline
    fn __radd__(self, other: Self.Coef) -> Self:
        '''Defines the reverse `+` plus operator. Returns other + self.'''
        return Self(self.re + other, self.im)

    @always_inline
    fn __rsub__(self, other: Self) -> Self:
        '''Defines the reverse `-` minus operator. Returns other - self.'''
        return other - self
    
    @always_inline
    fn __rsub__[__: None = None](self, other: Self.Lane) -> Self:
        '''Defines the reverse `-` minus operator. Returns other - self.'''
        return Self(other.re - self.re, other.im - self.im)

    @always_inline
    fn __rsub__(self, other: Self.Coef) -> Self:
        '''Defines the reverse `-` minus operator. Returns other - self.'''
        return Self(self.re - other, self.im)
    
    @always_inline
    fn __rmul__(self, other: Self) -> Self:
        '''Defines the reverse `*` product operator. Returns other * self.'''
        return other * self
    
    @always_inline
    fn __rmul__[__: None = None](self, other: Self.Lane) -> Self:
        '''Defines the reverse `*` product operator. Returns other * self.'''
        return Self(
            other.re * self.re - other.im * self.im, other.im * self.re + other.re * self.im
        )

    @always_inline
    fn __rmul__(self, other: Self.Coef) -> Self:
        '''Defines the reverse `*` product operator. Returns other * self.'''
        return Self(self.re * other, self.im * other)
    
    @always_inline
    fn __rtruediv__(self, other: Self) -> Self:
        '''Defines the reverse `/` divide operator. Returns other / self.'''
        return other / self
    
    @always_inline
    fn __rtruediv__[__: None = None](self, other: Self.Lane) -> Self:
        '''Defines the reverse `/` divide operator. Returns other / self.'''
        var denom: Self.Coef = 1 / self.squared_norm()
        return Self(
            (other.re * self.re + other.im * self.im) * denom, 
            (other.im * self.re - other.re * self.im) * denom,
        )

    @always_inline
    fn __rtruediv__(self, other: Self.Coef) -> Self:
        '''Defines the reverse `/` divide operator. Returns other / self.'''
        var denom: Self.Coef = 1 / self.squared_norm()
        return Self(self.re * other * denom, -self.im * other * denom)
    
    @always_inline
    fn __rfloordiv__(self, other: Self) -> Self:
        '''Defines the reverse `//` floor divide operator. Returns other // self.'''
        return other // self
    
    @always_inline
    fn __rfloordiv__[__: None = None](self, other: Self.Lane) -> Self:
        '''Defines the reverse `//` floor divide operator. Returns other // self.'''
        return round(other / self)
    
    @always_inline
    fn __rfloordiv__(self, other: Self.Coef) -> Self:
        '''Defines the reverse `//` floor divide operator. Returns other // self.'''
        return Self(other) // self
    
    @always_inline
    fn __rmod__(self, other: Self) -> Self:
        '''Defines the reverse `%` modulo operator. Returns other % self.'''
        return other % self

    @always_inline
    fn __rmod__[__: None = None](self, other: Self.Lane) -> Self:
        '''Defines the reverse `%` modulo operator. Returns other % self.'''
        return other - (other // self) * self

    @always_inline
    fn __rmod__(self, other: Self.Coef) -> Self:
        '''Defines the reverse `%` modulo operator. Returns other % self.'''
        return Self(other) % self

    @always_inline
    fn __rdivmod__(self, other: Self) -> Tuple[Self, Self]:
        '''Defines the reverse divmod operator. Returns the tuple (other // self, other % self).'''
        return other.__divmod__(self)
    
    @always_inline
    fn __rdivmod__[__: None = None](self, other: Self.Lane) -> Tuple[Self, Self]:
        '''Defines the reverse divmod operator. Returns the tuple (other // self, other % self).'''
        return (other // self, other % self)
    
    @always_inline
    fn __rdivmod__(self, other: Self.Coef) -> Tuple[Self, Self]:
        '''Defines the reverse divmod operator. Returns the tuple (other // self, other % self).'''
        return Self(other).__divmod__(self)
    
    @always_inline
    fn __rpow__(self, other: Self) -> Self:
        '''Defines the reverse `**` power operator. Returns other ** self.'''
        return other ** self

    @always_inline
    fn __rpow__[__: None = None](self, other: Self.Lane) -> Self:
        '''Defines the reverse `**` power operator. Returns other ** self.'''
        return Self(other.re, other.im) ** self
        
    @always_inline
    fn __rpow__(self, other: Self.Coef) -> Self:
        '''Defines the reverse `**` power operator. Returns other ** self.'''
        return Self(other) ** self

    # In-place arithmetic #############
    
    @always_inline
    fn __iadd__(mut self: Self, other: Self):
        '''Defines the `+=` in-place add operator. Computes self + other in-place.'''
        self.re += other.re
        self.im += other.im

    @always_inline
    fn __iadd__[__: None = None](mut self, other: Self.Coef):
        '''Defines the `+=` in-place add operator. Computes self + other in-place.'''
        self.re += other

    @always_inline
    fn __isub__(mut self, other: Self):
        '''Defines the `-=` in-place minus operator. Computes self - other in-place.'''
        self.re -= other.re
        self.im -= other.im

    @always_inline
    fn __isub__(mut self, other: Self.Coef):
        '''Defines the `-=` in-place minus operator. Computes self - other in-place.'''
        self.re -= other
    
    @always_inline
    fn __imul__(mut self, other: Self):
        '''Defines the `*=` in-place product operator. Computes self * other in-place.'''
        var new_re: Self.Coef = self.re * other.re - self.im * other.im
        self.im = self.re * other.im + self.im * other.re
        self.re = new_re

    @always_inline
    fn __imul__(mut self, other: Self.Coef):
        '''Defines the `*=` in-place product operator. Computes self * other in-place.'''
        self.re *= other
        self.im *= other
    
    @always_inline
    fn __imul__(mut self, other: Int):
        '''Defines the `*=` in-place product operator. Computes self * other in-place.'''
        self.re *= other
        self.im *= other

    @always_inline
    fn __itruediv__(mut self, other: Self):
        '''Defines the `/=` in-place divide operator. Computes self / other in-place.'''
        var denom: Self.Coef = 1 / other.squared_norm()
        var new_re: Self.Coef = (self.re * other.re + self.im * other.im) * denom
        self.im = (self.im * other.re - self.re * other.im) * denom
        self.re = new_re

    @always_inline
    fn __itruediv__(mut self, other: Self.Coef):
        '''Defines the `/=` in-place divide operator. Computes self / other in-place.'''
        self.re /= other
        self.im /= other
    
    @always_inline
    fn __ifloordiv__(mut self, other: Self):
        '''Defines the `//=` in-place floor divide operator. Computes self // other in-place.'''
        var quotient: Self = self / other
        self.re = round(quotient.re)
        self.im = round(quotient.im)
    
    @always_inline
    fn __ifloordiv__(mut self, other: Self.Coef):
        '''Defines the `//=` in-place floor divide operator. Computes self // other in-place.'''
        self.re = round(self.re / other)
        self.im = round(self.im / other)

    @always_inline
    fn __imod__(mut self, other: Self):
        '''Defines the `%=` in-place modulo operator. Computes self % other in-place.'''
        var mod: Self =  self - (self // other) * other
        self.re = mod.re
        self.im = mod.im
    
    @always_inline
    fn __imod__(mut self, other: Self.Coef):
        '''Defines the `%=` in-place modulo operator. Computes self % other in-place.'''
        var mod: Self =  self - (self // other) * other
        self.re = mod.re
        self.im = mod.im
    
    # Other ###########################

    @always_inline
    fn reciprocal(self) -> Self:
        '''Returns the reciprocal of a complex number.'''
        return self.conj() / self.squared_norm()

    @always_inline
    fn real(self) -> Self.Coef:
        '''Returns the real part of the complex number.'''
        return self.re

    @always_inline
    fn imag(self) -> Self.Coef:
        '''Returns the imaginary part of the complex number.'''
        return self.im
    
    @always_inline
    fn imaginary(self) -> Self.Coef:
        '''Alias of imag. Returns the real part of the complex number.'''
        return self.imag()
    
    @always_inline
    fn squared_norm(self) -> Self.Coef:
        '''Returns the squared norm a^2 + b^2 of a complex number a + bi.'''
        return self.re * self.re + self.im * self.im
    
    @always_inline
    fn norm(self) -> Self.Coef:
        '''Returns the absolute value of a complex number.'''
        return math.sqrt(self.squared_norm())
    
    @always_inline
    fn modulus(self) -> Self.Coef:
        '''Alias of norm. Returns the absolute value of a complex number.'''
        return self.norm()
    
    @always_inline
    fn arg(self) -> Self.Coef:
        '''Returns the argument atan2(b, a) of a complex number a + bi.'''
        return math.atan2(self.im, self.re)
    
    @always_inline
    fn phase(self) -> Self.Coef:
        '''Alias of arg. Returns the argument atan2(b, a) of a complex number a + bi.'''
        return self.arg()
    
    @always_inline
    fn conj(self) -> Self:
        '''Returns the conjugate a - bi of a complex number a + bi.'''
        return Self(self.re, -self.im)
    
    @always_inline
    fn conjugate(self) -> Self:
        '''Alias of conj. Returns the conjugate a - bi of a complex number a + bi.'''
        return self.conj()
    
    @always_inline
    fn exp(self) -> Self:
        '''Returns the exponential e^(a + bi) of a complex number a + bi.'''
        var r: Self.Coef = math.exp(self.re)
        return Self(r * math.cos(self.im), r * math.sin(self.im))
    
    @always_inline
    fn sqrt(self) -> Self:
        '''Returns the square root of a complex number using the principal branch.'''
        var r: Self.Coef = math.sqrt(self.norm())
        var theta: Self.Coef = self.arg() / 2
        return Self(r * math.cos(theta), r * math.sin(theta))

    @always_inline
    fn log(self) -> Self:
        '''Returns the natural logarithm of a complex number using the principal branch.'''
        return Self(math.log(self.norm()), self.arg())
    
    @always_inline
    fn log(self, base: Self.Coef) -> Self:
        '''Returns the logarithm in base `base` of a complex number using the principal branch.'''
        return self.log() / math.log(base)
    
    @always_inline
    fn sin(self) -> Self:
        '''Returns the sine of a complex number.'''
        return Self(math.sin(self.re) * math.cosh(self.im), math.cos(self.re) * math.sinh(self.im))

    @always_inline
    fn cos(self) -> Self:
        '''Returns the cosine of a complex number.'''
        return Self(
            math.cos(self.re) * math.cosh(self.im), -math.sin(self.re) * math.sinh(self.im)
        )
    
    @always_inline
    fn tan(self) -> Self:
        '''Returns the tangent of a complex number.'''
        return self.sin() / self.cos()
    
    @always_inline
    fn csc(self) -> Self:
        '''Returns the cosecant of a complex number.'''
        return self.sin().reciprocal()
    
    @always_inline
    fn sec(self) -> Self:
        '''Returns the secant of a complex number.'''
        return self.cos().reciprocal()
    
    @always_inline
    fn cot(self) -> Self:
        '''Returns the tangent of a complex number.'''
        return self.cos() / self.sin()
    
    @always_inline
    fn sinh(self) -> Self:
        '''Returns the hyperbolic sine of a complex number.'''
        return Self(
            math.sinh(self.re) * math.cos(self.im), math.cosh(self.re) * math.sin(self.im)
        )
    
    @always_inline
    fn cosh(self) -> Self:
        '''Returns the hyperbolic cosine of a complex number.'''
        return Self(
            math.cosh(self.re) * math.cos(self.im), math.sinh(self.re) * math.sin(self.im)
        )
    
    @always_inline
    fn tanh(self) -> Self:
        '''Returns the hyperbolic tangent of a complex number.'''
        return self.sinh() / self.cosh()
    
    @always_inline
    fn csch(self) -> Self:
        '''Returns the hyperbolic cosecant of a complex number.'''
        return self.sinh().reciprocal()
    
    @always_inline
    fn sech(self) -> Self:
        '''Returns the hyperbolic secant of a complex number.'''
        return self.cosh().reciprocal()
    
    @always_inline
    fn coth(self) -> Self:
        '''Returns the hyperbolic tangent of a complex number.'''
        return self.cosh() / self.sinh()
    
    @always_inline
    fn asin(self) -> Self:
        '''Returns the arcsin of a complex number, 
        enforcing the convention that the real part lies in [-pi/2, pi/2].
        '''
        var sqrt_term: Self = (Self(1) - self * self).sqrt()
        var result: Self = Self(0, -1) * (Self(-self.im, self.re) + sqrt_term).log()
        if result.re < -_halfpi:
            return Self(-result.re - math.pi, -result.im)
        if result.re > _halfpi:
            return Self(math.pi - result.re, -result.im)
        return result

    @always_inline
    fn acos(self) -> Self:
        '''Returns the arccos of a complex number, 
        enforcing the convention that the real part lies in [0, pi].
        '''
        var sqrt_term: Self = (self * self - Self(1)).sqrt()
        var result: Self = Self(0, -1) * (self + sqrt_term).log()
        if result.re < 0:
            result *= -1
        return result
    
    @always_inline
    fn atan(self) -> Self:
        '''Returns the arctan of a complex number, 
        enforcing the convention that the real part lies in [-pi/2, pi/2].
        '''
        var i: Self = Self(0, 1)
        return Self(0, -0.5) * ((i - self) / (i + self)).log()

    @always_inline
    fn asinh(self) -> Self:
        '''Returns the arcsinh of a complex number, 
        enforcing the convention that the imaginary part lies in [-pi/2, pi/2].
        '''
        var result: Self = (self + (self * self + Self(1)).sqrt()).log()
        if result.im < -_halfpi:
            return Self(-result.re, -result.im - math.pi)
        if result.im > _halfpi:
            return Self(-result.re, -result.im + math.pi)
        return result
    
    @always_inline
    fn acosh(self) -> Self:
        '''Returns the arccosh of a complex number, enforcing the convention that the 
        imaginary part lies in [-pi, pi] and the real part is non-negative.
        '''
        var result: Self = (self + (self * self - Self(1)).sqrt()).log()
        if result.re < 0:
            return Self(-result.re, -result.im)
        return result

    @always_inline
    fn atanh(self) -> Self:
        '''Returns the arctanh of a complex number, enforcing the convention that the
        imaginary part lies in [-pi/2, pi/2].
        '''
        var one: Self = Self(1)
        return Self(0.5) * ((one + self) / (one - self)).log()
        
    @always_inline
    fn floor(self) -> Self:
        '''Returns a complex number with the floor function applied to the real and
        imaginary parts.
        '''
        return Self(math.floor(self.re), math.floor(self.im))
    
    @always_inline
    fn ceil(self) -> Self:
        '''Returns a complex number with the ceiling function applied to the real and
        imaginary parts.
        '''
        return Self(math.ceil(self.re), math.ceil(self.im))
    
    @always_inline
    fn to_polar(self) -> Tuple[Self.Coef, Self.Coef]:
        '''Returns a tuple of the polar coordinates (r, θ) of a complex number r*e^(iθ).'''
        return (self.norm(), self.arg())
    
    @always_inline
    fn to_tuple(self) -> Tuple[Self.Coef, Self.Coef]:
        '''Returns the real and imaginary parts of a complex number as a tuple.'''
        return (self.re, self.im)
    
    @always_inline
    fn normalize(self) -> Self:
        '''Normalizes a complex number by dividing by its absolute value.'''
        return self / self.norm()
    
    @always_inline
    fn copy(self) -> Self:
        '''Return a copy of the complex number.'''
        return Self(self.re, self.im)
    
    @always_inline
    fn __hash__(self) -> UInt:
        '''Return a hash of the complex number. No idea if this hash function is any good...'''
        var real_prime: Int = 6355529899  # Random 10-digit prime
        var imag_prime: Int = 8422600973  # Random 10-digit prime
        
        var real_int: Int = Int(self.re * 1e8)
        var imag_int: Int = Int(self.im * 1e8)

        var real_hash: Int = real_int * real_prime
        var imag_hash: Int = imag_int * imag_prime
        return hash(real_hash ^ imag_hash)

    # Compare #########################
    
    @always_inline
    fn __eq__[__: None = None](self, other: Self) -> SIMD[DType.bool, Self.size]:
        '''Defines the `==` equality operator. Returns True if the real parts match and the 
        imaginary parts match.
        '''
        return self.re == other.re and self.im == other.im
    
    @always_inline
    fn __eq__(self, other: Self.Coef) -> SIMD[DType.bool, Self.size]:
        '''Defines the `==` equality operator. Returns True if the real parts match and the 
        imaginary parts match.
        '''
        return self.re == other and self.im == 0
    
    @always_inline
    fn __ne__[__: None = None](self, other: Self) -> SIMD[DType.bool, Self.size]:
        '''Defines the `!=` equality operator.
        Returns True if the real parts don't match or the imaginary parts don't match.
        '''
        return self.re != other.re or self.im != other.im
    
    @always_inline
    fn __ne__(self, other: Self.Coef) -> SIMD[DType.bool, Self.size]:
        '''Defines the `!=` equality operator.
        Returns True if the real parts don't match or the imaginary parts don't match.
        '''
        return self.re != other or self.im != 0
    
    @always_inline
    fn __lt__[__: None = None](self, other: Self) -> SIMD[DType.bool, Self.size]:
        '''Defines the `<` less than operator. 
        Returns True if the absolute value of self is less than the absolute value of other.
        '''
        return self.norm() < other.norm()
    
    @always_inline
    fn __lt__(self, other: Self.Coef) -> SIMD[DType.bool, Self.size]:
        '''Defines the `<` less than operator. Returns True if the absolute value of self is 
        less than other.
        '''
        return self.norm() < other
    
    @always_inline
    fn __le__[__: None = None](self, other: Self) -> SIMD[DType.bool, Self.size]:
        '''Defines the `<=` less than or equal operator. Returns True if the absolute value of 
        self is less than or equal to the absolute value of other.
        '''
        return self.norm() <= other.norm()
    
    @always_inline
    fn __le__(self, other: Self.Coef) -> SIMD[DType.bool, Self.size]:
        '''Defines the `<=` less than or equal operator. Returns True if the absolute value of
        self is less than or equal to other.
        '''
        return self.norm() <= other
    
    @always_inline
    fn __gt__[__: None = None](self, other: Self) -> SIMD[DType.bool, Self.size]:
        '''Defines the `<` greater than operator. Returns True if the absolute value of self is 
        greater than the absolute value of other.
        '''
        return self.norm() > other.norm()
    
    @always_inline
    fn __gt__(self, other: Self.Coef) -> SIMD[DType.bool, Self.size]:
        '''Defines the `>` greater than operator. Returns True if the absolute value of self is 
        greater than other.'''
        return self.norm() > other
    
    @always_inline
    fn __ge__[__: None = None](self, other: Self) -> SIMD[DType.bool, Self.size]:
        '''Defines the `>=` greater than or equal operator. Returns True if the absolute value of
        self is greater than or equal to the absolute value of other.
        '''
        return self.norm() >= other.norm()
    
    @always_inline
    fn __ge__(self, other: Self.Coef) -> SIMD[DType.bool, Self.size]:
        '''Defines the `>=` greater than or equal operator. Returns True if the absolute value of 
        self is greater than or equal to other.
        '''
        return self.norm() >= other
    
    @always_inline
    fn is_real(self) -> SIMD[DType.bool, Self.size]:
        '''Returns True if the imaginary component is zero.'''
        return self.im == 0
    
    @always_inline
    fn is_imag(self) -> SIMD[DType.bool, Self.size]:
        '''Returns True if the real component is zero.'''
        return self.re == 0
    
    @always_inline
    fn is_imaginary(self) -> SIMD[DType.bool, Self.size]:
        '''Alias of is_imag(). Returns True if the real component is zero.'''
        return self.is_imag()

    @always_inline
    fn is_zero(self) -> SIMD[DType.bool, Self.size]:
        '''Returns True if the real part and imaginaary part are zero.'''
        return self.re == 0 and self.im == 0
    
    @always_inline
    fn is_finite(self) -> SIMD[DType.bool, Self.size]:
        '''Returns True if the real and imaginary components are zero.'''
        return math.isfinite(self.re) and math.isfinite(self.im)
    
    @always_inline
    fn is_infinite(self) -> SIMD[DType.bool, Self.size]:
        '''Returns True if the real or imaginary part is infinite.'''
        return math.isinf(self.re) or math.isinf(self.im)
    
    @always_inline
    fn is_nan(self) -> SIMD[DType.bool, Self.size]:
        '''Returns True if the real or imaginary part is nan.'''
        return math.isnan(self.re) or math.isnan(self.im)
    
    @always_inline
    fn is_unit[tol: Scalar[Self.type] = DEFAULT_TOL](self) -> SIMD[DType.bool, Self.size]:
        '''Returns True if abs(self - 1) < tol.'''
        return abs(self - Self(1)) < tol
    
    @always_inline
    fn is_close[tol: Scalar[Self.type] = DEFAULT_TOL](
        self, other: Self
    ) -> SIMD[DType.bool, Self.size]:
        '''Returns True if abs(self - other) < tol.'''
        return abs(self - other) < tol
    
    @always_inline
    fn is_close[tol: Scalar[Self.type] = DEFAULT_TOL, __: None = None](
        self, other: Self.Lane
    ) -> Bool:
        '''Returns True if abs(self - other) < tol.'''
        return all(abs(self - other) < tol)
    
    @always_inline
    fn is_close[tol: Scalar[Self.type] = DEFAULT_TOL, __: None = None](
        self, other: Self.Coef
    ) -> SIMD[DType.bool, Self.size]:
        '''Returns True if abs(self - other) < tol.'''
        return abs(self - other) < tol
    
    @always_inline
    fn __len__(self) -> Int:
        return Self.size
    
    # Cast ############################
    
    @always_inline
    fn __bool__(self) -> SIMD[DType.bool, Self.size]:
        '''Returns true if teh complex number is not zero.'''
        return self.re != 0 or self.im != 0
    
    @always_inline
    fn __int__(self) -> Int:
        '''Casts the real part to an int.'''
        return Int(self.re)

    @always_inline
    fn __float__(self) -> Self.Coef:
        '''Casts the real part to a float.'''
        return self.re

    @always_inline
    fn __complex__(self) -> Self:
        '''Returns self (identity operation for complex numbers).'''
        return self
    
    @always_inline
    fn cast[target_dtype: DType](self) -> SIMD[target_dtype, Self.size]:
        '''Cast the real part to the specified DType.'''
        return self.re.cast[target_dtype]()

    # Item access #####################

    @always_inline
    fn __getitem__(self, idx: Int) -> Self.Lane:
        return Self.Lane(self.re[idx], self.im[idx])

    @always_inline
    fn __setitem__(mut self, idx: Int, item: Self.Lane):
        self.re[idx] = item.re
        self.im[idx] = item.im
    
    @always_inline
    fn slice[output_width: Int, /, *, offset: Int = 0](
        self
    ) -> ComplexSIMD[Self.type, output_width]:
        var re_slice: SIMD[Self.type, output_width] = self.re.slice[output_width, offset=offset]()
        var im_slice: SIMD[Self.type, output_width] = self.im.slice[output_width, offset=offset]()
        return ComplexSIMD[Self.type, output_width](re_slice, im_slice)
