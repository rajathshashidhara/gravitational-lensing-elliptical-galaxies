@everywhere using Cubature
using Iterators
using PyCall
pygui(:gtk)
using PyPlot

# defining global variables
@everywhere a = 1.0
@everywhere b = 1.0
@everywhere c = 1.0

acclevel = 1
st = 10

k = zeros(Float64, 2)

function cross2d(x::Array{Float64,1}, y::Array{Float64,1})
    x[1]*y[2] - x[2]*y[1]
end

function contains(aix::Array{Float64,1}, bix::Array{Float64,1}, cix::Array{Float64,1}, x::Array{Float64,1})
    da = x - aix
    db = x - bix
    dc = x - cix

    if cross2d(da,db)>=0 && cross2d(db, dc)>=0 && cross2d(dc,da)>=0
        return true
    elseif cross2d(da,db)<=0 && cross2d(db, dc)<=0 && cross2d(dc,da)<=0
        return true
    else
        return false
    end
end

@everywhere function sigma(x::Array{Float64,1})
    rho = 0.23873241463

    z = deepcopy(x)
    z[1] = z[1]/a
    z[2] = z[2]/b
    r = norm(z)

    if r<1.0
        2.0*c*rho*sqrt(1-r^2)
    else
        0.0
    end
end

@everywhere function dalpha(x::Array{Float64,1}, v::Array{Float64,1})
    m = ((k-x)/(norm(k-x))^2)*sigma(x)
    v[1] = m[1]
    v[2] = m[2]

    global k
    kd = k
    if isnan(v[1])
        k[1] += 1e-10
        dalpha(x,v)
    elseif isnan(isnan(v[2]))
        k[2] += 1e-10
        dalpha(x,v)
    end
    k = kd
end

@everywhere function quad(x::Array{Float64,1})
    global k
    k = x
    k - hcubature(2, dalpha, [-a;-b], [a;b]; reltol=1e-5,abstol=1e-5,maxevals=0)[1]
end

function testCubature()
    r = 0.25
    ϕ = zeros(Float64, 10)
    mag = zeros(Float64, 10)
    vect = []
    for θ in linspace(0.0, 2pi, 10)
        x = [r*cos(θ); r*sin(θ)]
        push!(vect, x)
    end

    rvec = pmap(quad, vect)
    for i=1:10
        mag[i] = norm(rvec[i])
        ϕ[i] = atan2(rvec[i][2], rvec[i][1])
    end
    plot(0.0:2pi:10, mag)
end

function imagesearch(pt::Array{Float64,1})
    bnd::Array{Array{Array{Float64,1},1},1} = Any[Any[[-5.0;5.0],[-5.0,5.0]]]
    fimages = cell(0)
    st = 24

    for k in 1:acclevel
        tbound = cell(0)
		images = cell(0)

        for bdry in bnd
            x = linspace(bdry[1][1], bdry[1][2], st)
            y = linspace(bdry[2][1], bdry[2][2], st)
            stepx = (bdry[1][2] - bdry[1][1])/(st-1)
			stepy = (bdry[2][2] - bdry[2][1])/(st-1)

            gd = []
            for ind in product(x,y)
                push!(gd, [ind...])
            end

            grid = pmap(quad,gd)

            for i in 1:(st-1)
				for j in 1:(st-1)
                    p = grid[i+(j-1)*st]
					q = grid[i+1+(j-1)*st]
					r = grid[i+1+j*st]
					s = grid[i+j*st]

                    if contains(p,q,r,pt) || contains(p,r,s,pt)
                        t = Any[[x[i]-0.1*stepx;x[i]+1*stepx],[y[j]-0.1*stepy, y[j]+1*stepy]]
                        push!(tbound, t)
                        push!(images, [x[i];y[j]])
                    end
                end
            end
        end
        bnd = tbound

        @printf "iteration: %d \t imagecount: %d \n" k length(bnd)
        @show bnd
        @printf "\n"

        fimages = images
    end
    @printf "\n"

    @printf "Images can be located at the following coordinates:\n"
    xim = []
    yim = []
    for im in fimages
        @printf "%+.8f %+.8f\n" im[1] im[2]
        push!(xim, im[1])
        push!(yim, im[2])
    end
    @printf "\n"

    fig = figure()
    scatter(xim, yim, color="red", marker="+")
    scatter([pt[1]], [pt[2]], color="blue", marker="o")
    savefig("images_$pt.svg")
end

imagesearch([0.3,0.1])
