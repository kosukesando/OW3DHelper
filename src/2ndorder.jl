# Calculate the second order corrections for eta & phi
using FFTW

g = 9.81
t = nothing
nx = oi.nx
ny = oi.ny
phase = deg2rad(oi.phase)
nkx = nothing
nky = nothing
kxi = nothing
kyi = nothing
depth = oi.depth
kxmatg = nothing
kymatg = nothing
kmatg = nothing
ωmatg = nothing
dirg = nothing
ampg_newwave_norm = nothing


C22 = zeros(2 * (nx - 1) + 1, 2 * (nx - 1) + 1)
V22 = zeros(2 * (nx - 1) + 1, 2 * (nx - 1) + 1)
C20 = zeros(2 * (nx - 1) + 1, 2 * (nx - 1) + 1)
V20 = zeros(2 * (nx - 1) + 1, 2 * (nx - 1) + 1)

ath = 10^(-6) # amplitude threshold for calculation of second order interactions

# iia=find(ampg_newwave_norm>ath)
# iia=ampg_newwave_norm>ath

#First calculate self interaction term.

maxes = 0

# for q00=1:1:length(iia)
for i in 1:nkx
    for j in 1:nky
        if ampg_newwave_norm[i, j] < ath
            continue
        end
        # [e2sx,iix2s]=min(abs(kxi-(kxvecg(ii)+kxvecg(ii))))
        # [e2sy,iiy2s]=min(abs(kyi-(kyvecg(ii)+kyvecg(ii))))
        e2sx = argmin(abs(kxi .- (kxmatg[i, j] + kxmatg[i, j])))
        iix2s = min(abs(kxi .- (kxmatg[i, j] + kxmatg[i, j])))
        e2sy = argmin(abs(kyi .- (kymatg[i, j] + kymatg[i, j])))
        iiy2s = min(abs(kyi .- (kymatg[i, j] + kymatg[i, j])))

        maxes = max([maxes e2sx e2sy])

        if maxes > (10^-12)
            "warning: misassigned component"
        end

        as2s = ((ampg_newwave_norm[i, j]^2) * kmatg[i, j] / (4 * tanh(kmatg[i, j] * depth))) * (2 + 3 / ((sinh(kmatg[i, j] * depth))^2))
        av2s = ((ampg_newwave_norm[i, j]^2) * (3 * ωmatg[i, j] / 8) * cosh(2 * kmatg[i, j] * depth) / ((sinh(kmatg[i, j] * depth))^4))
        p2s = (phase - t * ωmatg[i, j]) + (phase - t * ωmatg[i, j])
        c2s = complex(as2s * cos(p2s), as2s * sin(p2s))
        v2s = complex(av2s * cos(p2s - pi / 2), av2s * sin(p2s - pi / 2))

        C22[iiy2s, iix2s] = C22[iiy2s, iix2s] + c2s
        V22[iiy2s, iix2s] = V22[iiy2s, iix2s] + v2s
    end
end

#Now calculate the interaction between each component

println("Calculating second order cross-interactions...")

maxec = 0

# for q0 = 1:1:length(iia)
for ij1 in 1:nkx*nky-1
    i1 = 1 + (ij1 - 1) % nky
    j1 = 1 + (ij1 - 1) ÷ nky
    # ii = iia(q0)
    sinh1squared = (sinh(kmatg[i1, j1] * depth))^2
    tanh1 = tanh(kmatg[i1, j1] * depth)
    # for q1 = q0+1:1:length(iia) #Calculate the interaction only with the remaining terms
    for ij2 in ij1+1:nkx*nky
        i2 = 1 + (ij2 - 1) % nky
        j2 = 1 + (ij2 - 1) ÷ nky
        # jj = iia(q1)
        sinh2squared = (sinh(kmatg[i2, j2] * depth))^2
        tanh2 = tanh(kmatg[i2, j2] * depth)
        cosdt = cosd(dirg[i1, j1] - dirg[i2, j2]) # cos of the angle between the constituents
        pmodk = sqrt((kxmatg[i1, j1] + kxmatg[i2, j2])^2 + (kymatg[i1, j1] + kymatg[i2, j2])^2)  #mod(k1+k2)
        mmodk = sqrt((kxmatg[i1, j1] - kxmatg[i2, j2])^2 + (kymatg[i1, j1] - kymatg[i2, j2])^2)  #mod(k1-k2)
        Dp = (ωmatg[i1, j1] + ωmatg[i2, j2])^2 - g * pmodk * tanh(pmodk * depth) #Equation 21 in Dalzell
        Dm = (ωmatg[i1, j1] - ωmatg[i2, j2])^2 - g * mmodk * tanh(mmodk * depth) #Equation 22 in Dalzell
        tanhs = tanh1 * tanh2
        if tanhs < 0
            println("Alert! tanhs=$(tanhs)")
        end
        Ap1 = -(ωmatg[i1, j1] * ωmatg[i2, j2] * (ωmatg[i1, j1] + ωmatg[i2, j2]) / Dp) * (1 - cosdt / tanhs) #First term in eq 17 of Dalzell
        Ap2 = (1 / (2 * Dp)) * ((ωmatg[i1, j1]^3) / sinh1squared + (ωmatg[i2, j2]^3) / sinh2squared) #Second term in eq 17 of Dalzell
        Ap = Ap1 + Ap2
        Am1 = (ωmatg[i1, j1] * ωmatg[i2, j2] * (ωmatg[i1, j1] - ωmatg[i2, j2]) / Dm) * (1 + cosdt / tanhs) #First term in eq 18 of Dalzell
        Am2 = (1 / (2 * Dm)) * ((ωmatg[i1, j1]^3) / sinh1squared - (ωmatg[i2, j2]^3) / sinh2squared) #Second term in eq 18 of Dalzell
        Am = Am1 + Am2
        Bp0 = (ωmatg[i1, j1]^2 + ωmatg[i2, j2]^2) / (2 * g) #First term in 19 and 20 of Dalzell
        Bp1 = -((ωmatg[i1, j1] * ωmatg[i2, j2]) / (2 * g)) * (1 - cosdt / tanhs) * ((ωmatg[i1, j1] + ωmatg[i2, j2])^2 + g * pmodk * tanh(pmodk * depth)) / Dp #Second term in eq 19
        Bp2 = ((ωmatg[i1, j1] + ωmatg[i2, j2]) / (2 * g * Dp)) * (((ωmatg[i1, j1]^3) / sinh1squared) + ((ωmatg[i2, j2]^3) / sinh2squared)) #Third term in eq 19
        Bm1 = ((ωmatg[i1, j1] * ωmatg[i2, j2]) / (2 * g)) * (1 + cosdt / tanhs) * ((ωmatg[i1, j1] - ωmatg[i2, j2])^2 + g * mmodk * tanh(mmodk * depth)) / Dm #Second term in eq 20
        Bm2 = ((ωmatg[i1, j1] - ωmatg[i2, j2]) / (2 * g * Dm)) * (((ωmatg[i1, j1]^3) / sinh1squared) - ((ωmatg[i2, j2]^3) / sinh2squared)) #Third term in eq 19
        Bp = Bp0 + Bp1 + Bp2
        Bm = Bp0 + Bm1 + Bm2

        e2px = minimum(abs(kxi .- (kxmatg[i1, j1] + kxmatg[i2, j2])))
        e2py = minimum(abs(kyi .- (kymatg[i1, j1] + kymatg[i2, j2])))
        e2mx = minimum(abs(kxi .- (kxmatg[i2, j2] - kxmatg[i1, j1])))
        e2my = minimum(abs(kyi .- (kymatg[i2, j2] - kymatg[i1, j1])))
        iix2p = argmin(abs(kxi .- (kxmatg[i1, j1] + kxmatg[i2, j2])))
        iiy2p = argmin(abs(kyi .- (kymatg[i1, j1] + kymatg[i2, j2])))
        iix2m = argmin(abs(kxi .- (kxmatg[i2, j2] - kxmatg[i1, j1])))
        iiy2m = argmin(abs(kyi .- (kymatg[i2, j2] - kymatg[i1, j1])))

        maxec = maximum([maxec, e2px, e2py, e2mx, e2my])

        if maxec > (10^-12)
            println("warning: misassigned component")
        end

        a22 = ampg_newwave_norm[i1, j1] * ampg_newwave_norm[i2, j2]
        p22 = (phase - t * ωmatg[i1, j1]) + (phase - t * ωmatg[i2, j2])
        c22 = complex(a22 * Bp * cos(p22), a22 * Bp * sin(p22))
        v22 = complex(a22 * Ap * cos(p22 - pi / 2), a22 * Ap * sin(p22 - pi / 2))

        a20 = ampg_newwave_norm[i1, j1] * ampg_newwave_norm[i2, j2]
        p20 = (phase - t * ωmatg[i1, j1]) - (phase - t * ωmatg[i2, j2])
        c20 = complex(a20 * Bm * cos(p20), a20 * Bm * sin(p20))
        v20 = complex(a20 * Am * cos(p20 - pi / 2), a20 * Am * sin(p20 - pi / 2))

        C22[iiy2p, iix2p] = C22(iiy2p, iix2p) + c22
        V22[iiy2p, iix2p] = V22(iiy2p, iix2p) + v22
        C20[iiy2m, iix2m] = C20(iiy2m, iix2m) + c20
        V20[iiy2m, iix2m] = V20(iiy2m, iix2m) + v20
    end
end


# Construct Fourier Components corresponding to second order

F22 = 0.5 * prod(size(C22)) * ifftshift(ifftshift(C22, 1), 2)
S22 = circshift(circshift(ifftshift(ifftshift(ifft2(F22, "symmetric"), 1), 2), -1, 1), -1, 2)

G22 = 0.5 * prod(size(V22)) * ifftshift(ifftshift(V22, 1), 2)
P22 = circshift(circshift(ifftshift(ifftshift(ifft2(G22, "symmetric"), 1), 2), -1, 1), -1, 2)

F20 = 0.5 * prod(size(C20)) * ifftshift(ifftshift(C20, 1), 2)
S20 = circshift(circshift(ifftshift(ifftshift(ifft2(F20, "symmetric"), 1), 2), -1, 1), -1, 2)

G20 = 0.5 * prod(size(V20)) * ifftshift(ifftshift(V20, 1), 2)
P20 = fliplr(circshift(circshift(ifftshift(ifftshift(ifft2(G20, "symmetric"), 1), 2), -1, 1), -1, 2))
