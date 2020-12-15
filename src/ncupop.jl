mutable struct data_type
    supp
    basis
    coe
    obj
    tsupp
    ub
    sizes
end

function nctssos_first(f::Polynomial{false, T} where T<:Number, x::Vector{PolyVar{false}}; newton=true, reducebasis=true, TS="block", obj="eigen", merge=false, QUIET=false)
    n=length(x)
    mon=monomials(f)
    coe=coefficients(f)
    lm=length(mon)
    supp=[UInt16[] for i=1:lm]
    for i=1:lm
        ind=mon[i].z .>0
        vars=mon[i].vars[ind]
        exp=mon[i].z[ind]
        for j=1:length(vars)
            k=ncbfind(x, n, vars[j], rev=true)
            append!(supp[i], k*ones(UInt16, exp[j]))
        end
    end
    if obj=="trace"
        supp, coe=cyclic_canon(supp, coe)
    else
        supp, coe=sym_canon(supp, coe)
    end
    if newton==true
        if obj=="trace"
            d=Int(maximum(length.(supp))/2)
            basis=newton_cyclic(supp, n, d)
        else
            basis=newton_ncbasis(supp)
        end
    else
        basis=get_ncbasis(n, Int(maxdegree(f)/2))
    end
    tsupp=copy(supp)
    if obj=="trace"
        append!(tsupp, [_cyclic_canon([basis[i][end:-1:1]; basis[i]]) for i=1:length(basis)])
    else
        append!(tsupp, [[basis[i][end:-1:1]; basis[i]] for i=1:length(basis)])
    end
    sort!(tsupp)
    unique!(tsupp)
    blocks,cl,blocksize,ub,sizes,_=get_ncblocks(tsupp,basis,TS=TS,QUIET=QUIET,merge=merge,obj=obj)
    if reducebasis==true&&obj=="eigen"
        psupp=copy(supp)
        psupp=psupp[is_sym.(psupp)]
        push!(psupp, UInt16[])
        basis,flag=reducebasis!(psupp,basis,blocks,cl,blocksize)
        if flag==1
            tsupp=copy(supp)
            if obj=="trace"
                append!(tsupp, [_cyclic_canon([basis[i][end:-1:1]; basis[i]]) for i=1:length(basis)])
            else
                append!(tsupp, [[basis[i][end:-1:1]; basis[i]] for i=1:length(basis)])
            end
            sort!(tsupp)
            unique!(tsupp)
            blocks,cl,blocksize,ub,sizes,_=get_ncblocks(tsupp,basis,TS=TS,QUIET=QUIET,merge=merge,obj=obj)
        end
    end
    opt,tsupp=ncblockupop(supp,coe,basis,blocks,cl,blocksize,QUIET=QUIET,obj=obj)
    data=data_type(supp,basis,coe,obj,tsupp,ub,sizes)
    return opt,data
end

function nctssos_first(supp::Vector{Vector{UInt16}}, coe::Vector{Float64}, n::Int, d::Int; newton=true, reducebasis=true, TS="block", obj="eigen", merge=false, QUIET=false)
    if obj=="trace"
        supp, coe=cyclic_canon(supp, coe)
    else
        supp, coe=sym_canon(supp, coe)
    end
    if newton==true
        if obj=="trace"
            basis=newton_cyclic(supp, n, d)
        else
            basis=newton_ncbasis(supp)
        end
    else
        basis=get_ncbasis(n, d)
    end
    tsupp=copy(supp)
    if obj=="trace"
        append!(tsupp, [_cyclic_canon([basis[i][end:-1:1]; basis[i]]) for i=1:length(basis)])
    else
        append!(tsupp, [[basis[i][end:-1:1]; basis[i]] for i=1:length(basis)])
    end
    sort!(tsupp)
    unique!(tsupp)
    blocks,cl,blocksize,ub,sizes,_=get_ncblocks(tsupp,basis,TS=TS,QUIET=QUIET,merge=merge,obj=obj)
    if reducebasis==true&&obj=="eigen"
        psupp=copy(supp)
        psupp=psupp[is_sym.(psupp)]
        push!(psupp, UInt16[])
        basis,flag=reducebasis!(psupp,basis,blocks,cl,blocksize)
        if flag==1
            tsupp=copy(supp)
            if obj=="trace"
                append!(tsupp, [_cyclic_canon([basis[i][end:-1:1]; basis[i]]) for i=1:length(basis)])
            else
                append!(tsupp, [[basis[i][end:-1:1]; basis[i]] for i=1:length(basis)])
            end
            sort!(tsupp)
            unique!(tsupp)
            blocks,cl,blocksize,ub,sizes,_=get_ncblocks(tsupp,basis,TS=TS,QUIET=QUIET,merge=merge,obj=obj)
        end
    end
    opt,tsupp=ncblockupop(supp,coe,basis,blocks,cl,blocksize,QUIET=QUIET,obj=obj)
    data=data_type(supp,basis,coe,obj,tsupp,ub,sizes)
    return opt,data
end

function nctssos_higher!(data::data_type; TS="block", merge=false, QUIET=false)
    supp=data.supp
    basis=data.basis
    coe=data.coe
    obj=data.obj
    tsupp=data.tsupp
    ub=data.ub
    sizes=data.sizes
    blocks,cl,blocksize,ub,sizes,status=get_ncblocks(tsupp,basis,ub=ub,sizes=sizes,TS=TS,QUIET=QUIET,merge=merge,obj=obj)
    opt=nothing
    if status==1
        opt,tsupp=ncblockupop(supp,coe,basis,blocks,cl,blocksize,QUIET=QUIET,obj=obj)
    end
    data.tsupp=tsupp
    data.ub=ub
    data.sizes=sizes
    return opt,data
end

function cc(a::Vector{UInt16}, n::Int)
    ua=unique(a)
    ca=zeros(UInt8, n)
    for i=1:length(ua)
        ca[ua[i]]=count(x->isequal(ua[i], x), a)
    end
    return ca
end

function remove(csupp, dw, n)
    model=Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), true)
    t=@variable(model)
    alpha=@variable(model, [1:n])
    @constraint(model, [i=1:length(csupp)], alpha'*(csupp[i].-2*dw)<=t)
    @objective(model, Min, t)
    optimize!(model)
    if objective_value(model)>=0
        return true
    else
        return false
    end
end

function permutation(a)
    b=sparse(a)
    ua=convert(Vector{UInt16}, b.nzind)
    na=convert(Vector{UInt16}, b.nzval)
    return _permutation(ua, na)
end

function _permutation(ua, na)
    if !isempty(ua)
        perm=Vector{UInt16}[]
        for i=1:length(ua)
            nua=copy(ua)
            nna=copy(na)
            if na[i]==1
                deleteat!(nua, i)
                deleteat!(nna, i)
            else
                nna[i]-=1
            end
            temp=_permutation(nua, nna)
            push!.(temp, ua[i])
            append!(perm, temp)
        end
        return perm
    else
        return [UInt16[]]
    end
end

function newton_cyclic(supp, n, d)
    pbasis=get_basis(n,d)
    basis=[UInt16[]]
    csupp=cc.(supp, n)
    pushfirst!(csupp, zeros(UInt8, n))
    sort!(csupp)
    unique!(csupp)
    for i=2:size(pbasis,2)
        if remove(csupp, pbasis[:,i], n)
            append!(basis, permutation(pbasis[:,i]))
        end
    end
    sort!(basis)
    return basis
end

function get_basis(n,d)
    lb=binomial(n+d,d)
    basis=zeros(UInt8,n,lb)
    i=0
    t=1
    while i<d+1
        t+=1
        if basis[n,t-1]==i
           if i<d
              basis[1,t]=i+1
           end
           i+=1
        else
            j=findfirst(x->basis[x,t-1]!=0,1:n)
            basis[:,t]=basis[:,t-1]
            if j==1
               basis[1,t]-=1
               basis[2,t]+=1
            else
               basis[1,t]=basis[j,t]-1
               basis[j,t]=0
               basis[j+1,t]+=1
            end
        end
    end
    return basis
end

function _cyclic_canon(a::Vector{UInt16})
    if isempty(a)
        return a
    else
        return minimum([[a[i+1:length(a)]; a[1:i]] for i=0:length(a)-1])
    end
end

function cyclic_canon(supp, coe)
    nsupp=copy(supp)
    nsupp=_sym_canon.(nsupp)
    nsupp=_cyclic_canon.(nsupp)
    sort!(nsupp)
    unique!(nsupp)
    l=length(nsupp)
    ncoe=zeros(l)
    for i=1:length(supp)
        Locb=ncbfind(nsupp,l,_cyclic_canon(_sym_canon(supp[i])))
        ncoe[Locb]+=coe[i]
    end
    return nsupp, ncoe
end

function _sym_canon(a::Vector{UInt16})
    i=1
    while i<=Int(ceil((length(a)-1)/2))
        if a[i]<a[end+1-i]
            return a
        elseif a[i]>a[end+1-i]
            return reverse(a)
        else
            i+=1
        end
    end
    return a
end

function is_sym(a::Vector{UInt16})
    l=Int(ceil((length(a)-1)/2))
    return isequal(a[1:l], a[end:-1:end-l+1])
end

function sym_canon(supp, coe)
    nsupp=copy(supp)
    nsupp=_sym_canon.(nsupp)
    sort!(nsupp)
    unique!(nsupp)
    l=length(nsupp)
    ncoe=zeros(l)
    for i=1:length(supp)
        Locb=ncbfind(nsupp,l,_sym_canon(supp[i]))
        ncoe[Locb]+=coe[i]
    end
    return nsupp, ncoe
end

function get_ncbasis(n, d; ind=UInt16[i for i=1:n])
    basis=[UInt16[]]
    for i=1:d
        append!(basis, _get_ncbasis_deg(n, i, ind=ind))
    end
    return basis
end

function _get_ncbasis_deg(n, d; ind=UInt16[i for i=1:n])
    if d>0
        basis=Vector{UInt16}[]
        for i=1:n
            temp=_get_ncbasis_deg(n, d-1, ind=ind)
            push!.(temp, ind[i])
            append!(basis, temp)
        end
        return basis
    else
        return [UInt16[]]
    end
end

function newton_ncbasis(supp)
    nbasis=[UInt16[]]
    for bi in supp
        if iseven(length(bi))
            k=Int(length(bi)/2)
            w=bi[end-k+1:end]
            if isequal(w, bi[k:-1:1])
                for j=1:k
                    push!(nbasis, w[end-j+1:end])
                end
            end
        end
    end
    sort!(nbasis)
    unique!(nbasis)
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

function get_ncgraph(tsupp,basis;obj="eigen")
    lb=length(basis)
    G=SimpleGraph(lb)
    ltsupp=length(tsupp)
    for i = 1:lb, j = i+1:lb
        bi = [basis[i][end:-1:1]; basis[j]]
        bi=_sym_canon(bi)
        if obj=="trace"
            bi=_cyclic_canon(bi)
        end
        # if nx>0
        #     bi=comm(bi, nx)
        #     proj!(bi)
        # end
        if ncbfind(tsupp, ltsupp, bi)!=0
           add_edge!(G, i, j)
        end
    end
    return G
end

function get_ncblocks(tsupp,basis;ub=[],sizes=[],TS="block",obj="eigen",minimize=false,QUIET=true,merge=false)
    if TS==false
        blocksize=[length(basis)]
        blocks=[[i for i=1:length(basis)]]
        cl=1
    else
        G=get_ncgraph(tsupp,basis,obj=obj)
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
    tsupp=Vector{Vector{UInt16}}(undef, Int(sum(blocksize.^2+blocksize)/2))
    k=1
    for i=1:cl, j=1:blocksize[i], r=j:blocksize[i]
        @inbounds bi=[basis[blocks[i][j]][end:-1:1]; basis[blocks[i][r]]]
        @inbounds tsupp[k]=_sym_canon(bi)
        k+=1
    end
    if obj=="trace"
        tsupp=_cyclic_canon.(tsupp)
    end
    sort!(tsupp)
    unique!(tsupp)
    ltsupp=length(tsupp)
    model=Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons=[AffExpr(0) for i=1:ltsupp]
    for i=1:cl
        bs=blocksize[i]
        if bs==1
           @inbounds pos=@variable(model, lower_bound=0)
           @inbounds bi = [basis[blocks[i][1]][end:-1:1]; basis[blocks[i][1]]]
           if obj=="trace"
               bi=_cyclic_canon(bi)
           end
           Locb=ncbfind(tsupp,ltsupp,bi)
           @inbounds cons[Locb]+=pos
        else
           @inbounds pos=@variable(model, [1:bs, 1:bs], PSD)
           for j=1:blocksize[i], r=j:blocksize[i]
               @inbounds bi = [basis[blocks[i][j]][end:-1:1]; basis[blocks[i][r]]]
               bi=_sym_canon(bi)
               if obj=="trace"
                   bi=_cyclic_canon(bi)
               end
               Locb=ncbfind(tsupp,ltsupp,bi)
               if j==r
                   @inbounds cons[Locb]+=pos[j,r]
               else
                   @inbounds cons[Locb]+=2*pos[j,r]
               end
           end
        end
    end
    bc=zeros(ltsupp)
    for i=1:length(supp)
        Locb=ncbfind(tsupp,ltsupp,supp[i])
        if Locb==0
           @error "The monomial basis is not enough!"
           return nothing,nothing
        else
           bc[Locb]=coe[i]
        end
    end
    @variable(model, lower)
    cons[1]+=lower
    @constraint(model, con[i=1:ltsupp], cons[i]==bc[i])
    @objective(model, Max, lower)
    optimize!(model)
    status=termination_status(model)
    objv = objective_value(model)
    if status!=MOI.OPTIMAL
       println("termination status: $status")
       status=primal_status(model)
       println("solution status: $status")
    end
    println("optimum = $objv")
    return objv,tsupp
end

# function proj!(a::Vector{UInt16})
#     i=1
#     while i<=length(a)-1
#         if a[i]==a[i+1]
#             deleteat!(a, i)
#         else
#             i+=1
#         end
#     end
#     return a
# end
#
# function comm(a::Vector{UInt16}, nx)
#     ind1=a.<=nx
#     ind2=a.>nx
#     return [a[ind1]; a[ind2]]
# end
#
# function is_basis(a::Vector{UInt16}, nx)
#     for i=1:length(a)-1
#         if a[i]==a[i+1]
#             return false
#         end
#     end
#     return comm(a, nx)==a
# end
