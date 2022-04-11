using Singular

function add(l1 :: Vector{spoly{A}}, l2 :: Vector{spoly{A}}) where A
  length(l1) == length(l2) || error("vector dimensions not equal")
  ln = length(l1)
  R = parent(l1[1])
  carry = zero(R)
  l3 = [zero(R) for _ in 1:ln]
  for ix in 1:ln
    l3[ix] = carry+l1[ix]+l2[ix]
    #note that xy ∨ xz ∨ yz == xy + xz + yz
    carry = carry*l1[ix]+carry*l2[ix]+l1[ix]*l2[ix]
  end
  return l3
end

function sub(l1 :: Vector{spoly{A}}, l2 :: Vector{spoly{A}}) where A
  length(l1) == length(l2) || error("vector dimensions not equal")
  ln = length(l1)
  R = parent(l1[1])
  carry = one(R)
  l3 = [zero(R) for _ in 1:ln]
  for ix in 1:ln
    l3[ix] = carry+l1[ix]+l2[ix]+one(R)
    #note that xy ∨ xz ∨ yz == xy + xz + yz
    carry = carry*l1[ix]+carry*(l2[ix]+one(R))+l1[ix]*(l2[ix]+one(R))
  end
  return l3
end

function shr(u :: Vector{spoly{A}}; step :: Int = 1) where A
  length(u) > 0 || error("vector dimension equal to zero")
  R = parent(u[1])
  v = [zero(R) for _ in 1:length(u)]
  for ix in 1:length(u)-step
    v[ix] = u[ix+step]
  end
  return v
end

function and(l1 :: Vector{spoly{A}}, l2 :: Vector{spoly{A}}) where A
  length(l1) == length(l2) || error("vector dimensions not equal")
  ln = length(l1)
  R = parent(l1[1])
  l3 = [zero(R) for _ in 1:ln]
  for ix in 1:ln
    l3[ix] = l1[ix]*l2[ix]
  end
  return l3
end

#note: returns a vector of length `length(u)+step`
function shl_ext(u :: Vector{spoly{A}}; step :: Int = 1) where A
  length(u) > 0 || error("vector dimension equal to zero")
  R = parent(u[1])
  v = [zero(R) for _ in 1:(length(u)+step)]
  for ix in step+1:length(u)+step
    v[ix] = u[ix-step]
  end
  return v
end

#note: returns a vector of length `length(l1)+1`
function add_ext(l1 :: Vector{spoly{A}}, l2 :: Vector{spoly{A}}) where A
  length(l1) == length(l2) || error("vector dimensions not equal")
  ln = length(l1)
  R = parent(l1[1])
  carry = zero(R)
  l3 = [zero(R) for _ in 1:ln+1]
  for ix in 1:ln
    l3[ix] = carry+l1[ix]+l2[ix]
    #note that xy ∨ xz ∨ yz == xy + xz + yz
    carry = carry*l1[ix]+carry*l2[ix]+l1[ix]*l2[ix]
  end
  l3[ln+1] = carry
  return l3
end

function mapreduce(l::Vector{spoly{A}}, I::sideal{spoly{A}}) where A
  return map(x->reduce(x,I),l)
end

vars = ["b$ix" for ix in 0:15]
R, b = PolynomialRing(Fp(2), vars)
I = Ideal(R, [b[ix]^2-b[ix] for ix in 1:16])
I = std(I) #compute Gröbner basis of I
(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15) = b
reduce(b0^2+b1^3, I)

#0f = 00001111
#01 = 00000001
_5555 = [one(R),zero(R),one(R),zero(R),one(R),zero(R),one(R),zero(R),
         one(R),zero(R),one(R),zero(R),one(R),zero(R),one(R),zero(R)]
_3333 = [one(R),one(R),zero(R),zero(R),one(R),one(R),zero(R),zero(R),
         one(R),one(R),zero(R),zero(R),one(R),one(R),zero(R),zero(R)]
_0f0f = [one(R),one(R),one(R),one(R),zero(R),zero(R),zero(R),zero(R),
         one(R),one(R),one(R),one(R),zero(R),zero(R),zero(R),zero(R)]
g = [b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15]

di1 = deepcopy(g);
dx1 = shr(di1);
dx2 = and(dx1,_5555);

#di1 - dx2

di2 = sub(di1,dx2);
di3 = mapreduce(di2,I);
ax1 = deepcopy(di3);
di4 = shr(di3; step=2);
ax2 = and(ax1, _3333);
di5 = and(di4, _3333);
ax3 = mapreduce(ax2, I);
di6 = mapreduce(di5, I);
di7 = add(di6, ax3);
di8 = mapreduce(di7, I);
ax4 = deepcopy(di8);
ax5 = shr(ax4; step=4);
ax6 = add(ax5, di8);
ax7 = mapreduce(ax6, I);
ax8 = and(ax7, _0f0f);
ax9 = mapreduce(ax8, I);
#imul   ax,ax,0x0101
#0x0101 = 0b0000000100000001
ax10 = [ax9[1],ax9[2],ax9[3],ax9[4],ax9[5],ax9[6],ax9[7],ax9[8],
         ax9[9],ax9[10],ax9[11],ax9[12],ax9[13],ax9[14],ax9[15],ax9[16],
         zero(R),zero(R),zero(R),zero(R),zero(R),zero(R),zero(R),zero(R)];
ax11 = shl_ext(ax9; step=8);
ax12 = add_ext(ax10,ax11);
ax13 = mapreduce(ax12, I)
#shr ax, 0x8
ax14 = [ax13[9],ax13[10],ax13[11],ax13[12],
        ax13[13],ax13[14],ax13[15],ax13[16],
        zero(R),zero(R),zero(R),zero(R),
        zero(R),zero(R),zero(R),zero(R)]

function equiv(u::BitVector)
  n = sum(a)
  v = [0 for _ in 1:32]
  v[1] = n % 2
  ix = 1
  while 2^ix <= n
    v[ix+1] = binomial(n,2^ix) % 2
    ix += 1
  end
  return v
end
