using Singular
using Singular: mul!, addeq!, zero!

function shreq!(v :: Vector{spoly{A}}, u :: Vector{spoly{A}}; step :: Int = 1) where A
  (length(v) > 0 && length(u) > 0) || error("vector dimension equal to zero")
  R = parent(v[1])
  for ix in 1:length(u)-step
    v[ix] = u[ix+step]
  end
end


function _addeq!(l0 :: Vector{spoly{A}},
                 l1 :: Vector{spoly{A}},
                 l2 :: Vector{spoly{A}},
                 I  :: sideal{spoly{A}}) where A
  length(l0) == length(l1) == length(l2) || error("vector dimensions not equal")
  ln = length(l0)
  R = parent(l0[1])
  temp1 = zero(R)
  temp2 = zero(R)
  temp3 = zero(R)
  l0_ix_unred = zero(R)
  carry_unred = zero(R)
  carry = zero(R)
  for ix in 1:ln
    #w/o inlining -- l0[ix] = carry+l1[ix]+l2[ix]
    zero!(l0_ix_unred)
    addeq!(l0_ix_unred, l1[ix])
    addeq!(l0_ix_unred, l2[ix])
    addeq!(l0_ix_unred, carry)
    l0[ix] = reduce(l0_ix_unred, I)
    #note that xy ∨ xz ∨ yz == xy + xz + yz
    #
    #w/o inlining -- carry = carry*l1[ix]+carry*l2[ix]+l1[ix]*l2[ix]
    mul!(temp1, l1[ix], l2[ix])
    mul!(temp2, carry, l2[ix])
    mul!(temp3, carry, l1[ix])
    zero!(carry_unred)
    addeq!(carry_unred, temp1)
    addeq!(carry_unred, temp2)
    addeq!(carry_unred, temp3)
    zero!(carry)
    carry = reduce(carry_unred, I)
  end
end

function _addeq!(l0 :: Vector{spoly{A}},
                 l1 :: Vector{spoly{A}},
                 l2 :: Vector{spoly{A}}) where A
  length(l0) == length(l1) == length(l2) || error("vector dimensions not equal")
  ln = length(l0)
  R = parent(l0[1])
  temp1 = zero(R)
  temp2 = zero(R)
  temp3 = zero(R)
  l0_ix_unred = zero(R)
  carry_unred = zero(R)
  carry = zero(R)
  for ix in 1:ln
    #w/o inlining -- l0[ix] = carry+l1[ix]+l2[ix]
    addeq!(l0[ix], l1[ix])
    addeq!(l0[ix], l2[ix])
    addeq!(l0[ix], carry)
    #note that xy ∨ xz ∨ yz == xy + xz + yz
    #
    #w/o inlining -- carry = carry*l1[ix]+carry*l2[ix]+l1[ix]*l2[ix]
    mul!(temp1, l1[ix], l2[ix])
    mul!(temp2, carry, l2[ix])
    mul!(temp3, carry, l1[ix])
    zero!(carry)
    addeq!(carry, temp1)
    addeq!(carry, temp2)
    addeq!(carry, temp3)
  end
end


function subeq!(l0 :: Vector{spoly{A}}, 
                l1 :: Vector{spoly{A}}, 
                l2 :: Vector{spoly{A}}, 
                I  :: sideal{spoly{A}}) where A
  length(l0) == length(l1) == length(l2) || error("vector dimensions not equal")
  ln = length(l0)
  R = parent(l0[1])
  temp1 = zero(R)
  temp2 = zero(R)
  temp3 = zero(R)
  l0_ix_unred = zero(R)
  carry_unred = one(R)
  carry = one(R)
  for ix in 1:ln
    #w/o inlining -- l0[ix] = carry+l1[ix]+l2[ix]+one(R)
    zero!(l0_ix_unred)
    addeq!(l0_ix_unred, one(R))
    addeq!(l0_ix_unred, l2[ix])
    addeq!(l0_ix_unred, l1[ix])
    addeq!(l0_ix_unred, carry)
    l0[ix] = reduce(l0_ix_unred, I)
    #note that xy ∨ xz ∨ yz == xy + xz + yz
    #
    #w/o inlining -- carry = carry*l1[ix]+carry*l2[ix]+l1[ix]*l2[ix]+l1[ix]+carry
    mul!(temp1, l1[ix], l2[ix])
    mul!(temp2, carry, l2[ix])
    mul!(temp3, carry, l1[ix])
    zero!(carry_unred)
    addeq!(carry_unred, carry)
    addeq!(carry_unred, l1[ix])
    addeq!(carry_unred, temp1)
    addeq!(carry_unred, temp2)
    addeq!(carry_unred, temp3)
    zero!(carry)
    carry = reduce(carry_unred, I)
  end
end

function andeq!(l0 :: Vector{spoly{A}},
                l1 :: Vector{spoly{A}},
                l2 :: Vector{spoly{A}},
                I  :: sideal{spoly{A}}) where A
  length(l0) == length(l1) == length(l2) || error("vector dimensions not equal")
  ln = length(l0)
  R = parent(l0[1])
  l0_ix_unred = zero(R)
  for ix in 1:ln
    zero!(l0_ix_unred)
    mul!(l0_ix_unred,l1[ix],l2[ix])
    l0[ix] = reduce(l0_ix_unred, I)
  end
end

vars = ["b$ix" for ix in 0:31]
R, b = PolynomialRing(Fp(2), vars)
I = Ideal(R, [b[ix]^2-b[ix] for ix in 1:32])
I = std(I) #compute Gröbner basis of I
(b0,b1,b2,b3,b4,b5,b6,b7,
 b8,b9,b10,b11,b12,b13,b14,b15,
 b16,b17,b18,b19,b20,b21,b22,b23,
 b24,b25,b26,b27,b28,b29,b30,b31) = b
reduce(b0^2+b1^3, I)


_55555555 = [one(R),zero(R),one(R),zero(R),one(R),zero(R),one(R),zero(R),
             one(R),zero(R),one(R),zero(R),one(R),zero(R),one(R),zero(R),
             one(R),zero(R),one(R),zero(R),one(R),zero(R),one(R),zero(R),
             one(R),zero(R),one(R),zero(R),one(R),zero(R),one(R),zero(R)]

_33333333 = [one(R),one(R),zero(R),zero(R),one(R),one(R),zero(R),zero(R),
             one(R),one(R),zero(R),zero(R),one(R),one(R),zero(R),zero(R),
             one(R),one(R),zero(R),zero(R),one(R),one(R),zero(R),zero(R),
             one(R),one(R),zero(R),zero(R),one(R),one(R),zero(R),zero(R)]

_0f0f0f0f = [one(R),one(R),one(R),one(R),zero(R),zero(R),zero(R),zero(R),
             one(R),one(R),one(R),one(R),zero(R),zero(R),zero(R),zero(R),
             one(R),one(R),one(R),one(R),zero(R),zero(R),zero(R),zero(R),
             one(R),one(R),one(R),one(R),zero(R),zero(R),zero(R),zero(R)]

_01010101 = [one(R),zero(R),zero(R),zero(R),zero(R),zero(R),zero(R),zero(R),
             one(R),zero(R),zero(R),zero(R),zero(R),zero(R),zero(R),zero(R),
             one(R),zero(R),zero(R),zero(R),zero(R),zero(R),zero(R),zero(R),
             one(R),zero(R),zero(R),zero(R),zero(R),zero(R),zero(R),zero(R)]
#mov    edx,edi
#shr    edx,1
#and    edx,0x55555555
#sub    edi,edx
#mov    eax,edi
#shr    edi,0x2
#and    eax,0x33333333
#and    edi,0x33333333
#add    edi,eax
#mov    eax,edi
#shr    eax,0x4
#add    eax,edi
#and    eax,0x0f0f0f0f

edi1 = [b0,b1,b2,b3,b4,b5,b6,b7,
        b8,b9,b10,b11,b12,b13,b14,b15,
        b16,b17,b18,b19,b20,b21,b22,b23,
        b24,b25,b26,b27,b28,b29,b30,b31]
edx1 = [zero(R) for _ in 1:32]
shreq!(edx1, edi1);
edx2 = [zero(R) for _ in 1:32]
andeq!(edx2, edx1,_55555555, I);
edi2 = [zero(R) for _ in 1:32]
subeq!(edi2, edi1, edx2, I);
eax1 = deepcopy(edi2);
edi3 = [zero(R) for _ in 1:32]
shreq!(edi3, edi2; step=2);
eax2 = [zero(R) for _ in 1:32]
andeq!(eax2, eax1, _33333333, I)
edi4 = [zero(R) for _ in 1:32]
andeq!(edi4, edi3, _33333333, I)
edi5 = [zero(R) for _ in 1:32]
_addeq!(edi5, edi4, eax2, I)
eax3 = deepcopy(edi5)
eax4 = [zero(R) for _ in 1:32]
shreq!(eax4, eax3, step=4);
eax5 = [zero(R) for _ in 1:32]
_addeq!(eax5, eax4, edi5, I);
eax6 = [zero(R) for _ in 1:32]
andeq!(eax6, eax5, _0f0f0f0f, I);


#time for imul

vars1 = ["b0_7_1","b0_7_2","b0_7_4","b0_7_8",
         "b8_15_1","b8_15_2","b8_15_4","b8_15_8",
         "b16_23_1","b16_23_2","b16_23_4","b16_23_8",
         "b24_31_1","b24_31_2","b24_31_4","b24_31_8"]
R1, x = PolynomialRing(Fp(2), vars1)
I1 = Ideal(R1, [x[ix]^2-x[ix] for ix in 1:16])
I1 = std(I1)
(x0,x1,x2,x3,x4,x5,x6,x7,
 x8,x9,x10,x11,x12,x13,x14,x15) = x

(b0_7_1,b0_7_2,b0_7_4,b0_7_8,
 b8_15_1,b8_15_2,b8_15_4,b8_15_8,
 b16_23_1,b16_23_2,b16_23_4,b16_23_8,
 b24_31_1,b24_31_2,b24_31_4,b24_31_8) = x

v1 = [x0,x1,x2,x3,zero(R1),zero(R1),zero(R1),zero(R1),
      x4,x5,x6,x7,zero(R1),zero(R1),zero(R1),zero(R1),
      x8,x9,x10,x11,zero(R1),zero(R1),zero(R1),zero(R1),
      x12,x13,x14,x15,zero(R1),zero(R1),zero(R1),zero(R1)]
v2 = [zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),
      x0,x1,x2,x3,zero(R1),zero(R1),zero(R1),zero(R1),
      x4,x5,x6,x7,zero(R1),zero(R1),zero(R1),zero(R1),
      x8,x9,x10,x11,zero(R1),zero(R1),zero(R1),zero(R1)]
v3 = [zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),
      zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),
      x0,x1,x2,x3,zero(R1),zero(R1),zero(R1),zero(R1),
      x4,x5,x6,x7,zero(R1),zero(R1),zero(R1),zero(R1)]
v4 = [zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),
      zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),
      zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),zero(R1),
      x0,x1,x2,x3,zero(R1),zero(R1),zero(R1),zero(R1)]


#imul   eax,eax,0x01010101
#shr    eax,0x18
#ret

v5 = [zero(R1) for _ in 1:32]
_addeq!(v5,v1,v2,I1)
v6 = [zero(R1) for _ in 1:32]
_addeq!(v6,v5,v3,I1)
v7 = [zero(R1) for _ in 1:32]
_addeq!(v7,v6,v4,I1)
v8 = [zero(R1) for _ in 1:32]
shreq!(v8,v7;step=24)

@assert(all(v8[ix]==zero(R1) for ix in 7:32))

#note that
#
#b0_7_1 = all 1-factor products                 with factors in b0,b1,b2,...,b7
#b0_7_2 = all 2-factor products
#b0_7_1*b0_7_2 = all 3-factor products
#b0_7_4 = all 4-factor products
#b0_7_1*b0_7_4 = all 5-factor products
#b0_7_2*b0_7_4 = all 6-factor products
#b0_7_1*b0_7_2*b0_7_4 = all 7-factor products
#b0_7_8 = all 8-factor products

#hence v8[2] is sum of all 2-factor products    with factors in b0,b1,b2,...,b31
#      v8[3] is sum of all 4-factor products
#      v8[4] is sum of all 8-factor products
#      v8[5] is sum of all 16-factor products
#      v8[6] is sum of all 32-factor products
