mutable struct data_type
    supp
    basis
    coe
    obj
    supp1
    ub
    sizes
end

function ncblockupop_first(f, x; newton=true, reducebasis=true, TS="block", obj="eigen", merge=false, QUIET=false)
    n=length(x)
    if obj=="trace"
        f=cyclic_canon(f, x)
    end
    mon=monomials(f)
    coe=coefficients(f)
    lm=length(mon)
    supp=[UInt16[] for i=1:lm]
    for i=1:lm
        vars=variables(mon[i])
        exp=exponents(mon[i])
        ind=[exp[k]!=0 for k=1:length(exp)]
        vars=vars[ind]
        exp=exp[ind]
        for j=1:length(vars)
            k=ncbfind(x, n, vars[j], rev=true)
            append!(supp[i], k*ones(UInt16, exp[j]))
        end
    end
    d=Int(maxdegree(f)/2)
    if newton==true
        basis=newton_ncbasis(n, d, supp)
    else
        basis=get_ncbasis(n, d)
    end
    blocks,cl,blocksize,ub,sizes,_=get_ncblocks(supp,basis,TS=TS,QUIET=QUIET,merge=merge,obj=obj)
    if reducebasis==true
        tsupp=copy(supp)
        push!(tsupp, UInt16[])
        basis,flag=reducebasis!(tsupp,basis,blocks,cl,blocksize)
        if flag==1
            blocks,cl,blocksize,ub,sizes,_=get_ncblocks(supp,basis,TS=TS,QUIET=QUIET,merge=merge,obj=obj)
        end
    end
    opt,supp1=ncblockupop(supp,coe,basis,blocks,cl,blocksize,QUIET=QUIET,obj=obj)
    data=data_type(supp,basis,obj,coe,supp1,ub,sizes)
    return opt,data
end

function ncblockupop_higher!(data;TS="block",merge=false,QUIET=false)
    supp=data.supp
    basis=data.basis
    coe=data.coe
    obj=data.obj
    supp1=data.supp1
    ub=data.ub
    sizes=data.sizes
    blocks,cl,blocksize,ub,sizes,status=get_ncblocks(supp,basis,ub=ub,sizes=sizes,TS=TS,QUIET=QUIET,merge=merge,obj=obj)
    opt=nothing
    if status==1
        opt,supp1=ncblockupop(supp,coe,basis,blocks,cl,blocksize,QUIET=QUIET,obj=obj)
    end
    data.supp1=supp1
    data.ub=ub
    data.sizes=sizes
    return opt,data
end

function cyclic_canon(f, x)
    mon=monomials(f)
    coe=coefficients(f)
    return sum([coe[i]*prod(x[_cyclic_canon(mon[i], x)]) for i=1:length(mon)])
end

function _cyclic_canon(w, x)
    if isconstant(w)
        return UInt16[]
    end
    vars=variables(w)
    exp=exponents(w)
    ind=[exp[i]!=0 for i=1:length(exp)]
    vars=vars[ind]
    exp=exp[ind]
    a=UInt16[]
    for j=1:length(vars)
        k=ncbfind(x, n, vars[j], rev=true)
        append!(a, k*ones(UInt16, exp[j]))
    end
    return minimum([[a[i+1:length(a)]; a[1:i]] for i=0:length(a)-1])
end

function _cyclic_canon(a)
    if isempty(a)
        return UInt16[]
    else
        return minimum([[a[i+1:length(a)]; a[1:i]] for i=0:length(a)-1])
    end
end

function symmetrize(f)
    mon=monomials(f)
    coe=coefficients(f)
    g=sum([coe[i]*prod(variables(mon[i])[end:-1:1].^exponents(mon[i])[end:-1:1]) for i=1:length(mon)])
    return (f+g)/2
end

function iscyclic(a, b)
    n=length(a)
    if n!=length(b)
        return false
    else
        for i=0:n-1
            if isequal([a[i+1:n]; a[1:i]], b)
                return true
            end
        end
        return false
    end
end

function cyclic_class(A)
    n=length(A)
    G=SimpleGraph(n)
    for i = 1:lb, j = i+1:lb
        if iscyclic(A[i], A[j])
           add_edge!(G, i, j)
        end
    end
    blocks=connected_components(G)
    blocksize=length.(blocks)
    cl=length(blocksize)
    return blocks,cl,blocksize
end

function get_ncbasis(n,d)
    basis=[UInt16[]]
    for i=1:d
        append!(basis, _get_ncbasis_deg(n,i))
    end
    return basis
end

function _get_ncbasis_deg(n,d)
    if d>0
        basis=Vector{UInt16}[]
        for i=1:n
            temp=_get_ncbasis_deg(n,d-1)
            push!.(temp, i)
            append!(basis, temp)
        end
        return basis
    else
        return [UInt16[]]
    end
end

function newton_ncbasis(n, d, supp)
    tsupp=copy(supp)
    sort!(tsupp)
    ltsupp=length(tsupp)
    basis=get_ncbasis(n, d)
    lb=length(basis)
    nbasis=[UInt16[]]
    for i=2:lb
        bi=[basis[i][end:-1:1]; basis[i]]
        if ncbfind(tsupp, ltsupp, bi)!=0
            for j=1:length(basis[i])
                temp=basis[i][end-j+1:end]
                push!(nbasis, temp)
            end
        end
    end
    sort!(nbasis)
    return nbasis
end

function ncbfind(A, l, a; rev=false)
    if l==0
        return 0
    end
    low=1
    high=l
    while low<=high
        mid=Int(ceil(1/2*(low+high)))
        if isequal(A[mid], a)
           return mid
        elseif isless(A[mid], a)
            if rev==false
                low=mid+1
            else
                high=mid-1
            end
        else
            if rev==false
                high=mid-1
            else
                low=mid+1
            end
        end
    end
    return 0
end

function get_ncgraph(supp,basis;obj="eigen")
    lb=length(basis)
    G=SimpleGraph(lb)
    tsupp=copy(supp)
    if obj=="trace"
        append!(tsupp, [_cyclic_canon([basis[i][end:-1:1]; basis[i]]) for i=1:lb])
    else
        append!(tsupp, [[basis[i][end:-1:1]; basis[i]] for i=1:lb])
    end
    sort!(tsupp)
    unique!(tsupp)
    ltsupp=length(tsupp)
    for i = 1:lb, j = i+1:lb
        bi = [basis[i][end:-1:1]; basis[j]]
        if obj=="trace"
            bi=_cyclic_canon(bi)
        end
        if ncbfind(tsupp, ltsupp, bi)!=0
           add_edge!(G, i, j)
        end
    end
    return G, tsupp
end

function get_ncblocks(supp,basis;ub=[],sizes=[],TS="block",obj="eigen",minimize=false,QUIET=true,merge=false)
    if TS==false
        blocksize=[length(basis)]
        blocks=[[i for i=1:length(basis)]]
        cl=1
    else
        G,_=get_ncgraph(supp,basis,obj=obj)
        if TS=="block"
            blocks=connected_components(G)
            blocksize=length.(blocks)
            cl=length(blocksize)
        else
            blocks,cl,blocksize=chordal_cliques!(G, method=TS, minimize=minimize)
            if merge==true
                blocks,cl,blocksize=clique_merge!(blocks,cl,QUIET=true)
            end
        end
    end
    nub=unique(blocksize)
    nsizes=[sum(blocksize.== i) for i in nub]
    if isempty(ub)||nub!=ub||nsizes!=sizes
        status=1
        if QUIET==false
            println("------------------------------------------------------")
            println("The sizes of blocks:\n$nub\n$nsizes")
            println("------------------------------------------------------")
        end
    else
        status=0
        if QUIET==false
            println("No higher NCTSSOS hierarchy!")
        end
    end
    return blocks,cl,blocksize,nub,nsizes,status
end

function ncblockupop(supp,coe,basis,blocks,cl,blocksize;QUIET=true,obj="eigen")
    supp1=Vector{Vector{UInt16}}(undef, sum(blocksize.^2))
    k=1
    for i=1:cl, j=1:blocksize[i], r=1:blocksize[i]
        @inbounds bi=[basis[blocks[i][j]][end:-1:1]; basis[blocks[i][r]]]
        @inbounds supp1[k]=bi
        k+=1
    end
    if obj=="trace"
        supp1=_cyclic_canon.(supp1)
    end
    sort!(supp1)
    unique!(supp1)
    lsupp1=length(supp1)
    model=Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons=[AffExpr(0) for i=1:lsupp1]
    for i=1:cl
        bs=blocksize[i]
        if bs==1
           @inbounds pos=@variable(model, lower_bound=0)
           @inbounds bi = [basis[blocks[i][1]][end:-1:1]; basis[blocks[i][1]]]
           if obj=="trace"
               bi=_cyclic_canon(bi)
           end
           Locb=ncbfind(supp1,lsupp1,bi)
           @inbounds cons[Locb]+=pos
        else
           @inbounds pos=@variable(model, [1:bs, 1:bs], PSD)
           for j=1:blocksize[i], r=1:blocksize[i]
               @inbounds bi = [basis[blocks[i][j]][end:-1:1]; basis[blocks[i][r]]]
               if obj=="trace"
                   bi=_cyclic_canon(bi)
               end
               Locb=ncbfind(supp1,lsupp1,bi)
               @inbounds cons[Locb]+=pos[j,r]
           end
        end
    end
    bc=zeros(lsupp1)
    lsupp=length(supp)
    for i=1:lsupp
        Locb=ncbfind(supp1,lsupp1,supp[i])
        if Locb==0
           @error "The monomial basis is not enough!"
           return nothing,nothing
        else
           bc[Locb]=coe[i]
        end
    end
    @variable(model, lower)
    cons[1]+=lower
    @constraint(model, con[i=1:lsupp1], cons[i]==bc[i])
    @objective(model, Max, lower)
    optimize!(model)
    status=termination_status(model)
    if status == MOI.OPTIMAL
       objv = objective_value(model)
       println("optimum = $objv")
    else
       objv = objective_value(model)
       println("termination status: $status")
       sstatus=primal_status(model)
       println("solution status: $sstatus")
       println("optimum = $objv")
    end
    return objv,supp1
end
