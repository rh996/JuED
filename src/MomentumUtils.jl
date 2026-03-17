module MomentumUtilsMod

export momentum_add_1d, momentum_sub_1d, momentum_add_2d, momentum_sub_2d

@inline function momentum_add_1d(k1, k2, systemsize)
    return mod(k1 + k2, systemsize)
end

@inline function momentum_sub_1d(k1, k2, systemsize)
    return mod(k1 - k2, systemsize)
end

@inline function momentum_add_2d(k1, k2, Nkx, Nky)
    k1x = mod(k1, Nkx)
    k1y = fld(k1, Nkx)
    k2x = mod(k2, Nkx)
    k2y = fld(k2, Nkx)
    return mod(k1x + k2x, Nkx) + mod(k1y + k2y, Nky) * Nkx
end

@inline function momentum_sub_2d(k1, k2, Nkx, Nky)
    k1x = mod(k1, Nkx)
    k1y = fld(k1, Nkx)
    k2x = mod(k2, Nkx)
    k2y = fld(k2, Nkx)
    return mod(k1x - k2x, Nkx) + mod(k1y - k2y, Nky) * Nkx
end

end
