mutable struct cdata_type
    m
    numeq
    supp
    coe
    obj
    lt
    fbasis
    gbasis
    tsupp
    ub
    sizes
    gblocks
    gcl
    gblocksize
end

function nctssos_first(pop::Vector{Polynomial{false, T}} where T<:Number,x::Vector{PolyVar{false}},d::Int;numeq=0,reducebasis=false,TS="block",obj="eigen",merge=false,QUIET=false)
    n=length(x)
    m=length(pop)-1
    coe=Vector{Vector{Float64}}(undef, m+1)
    supp=Vector{Vector{Vector{UInt16}}}(undef, m+1)
    lt=Vector{UInt16}(undef, m+1)
    for k=1:m+1
        mon=monomials(pop[k])
        coe[k]=coefficients(pop[k])
        lt[k]=length(mon)
        supp[k]=[UInt16[] for i=1:lt[k]]
        for i=1:lt[k]
            ind=mon[i].z .>0
            vars=mon[i].vars[ind]
            exp=mon[i].z[ind]
            for j=1:length(vars)
                l=ncbfind(x, n, vars[j], rev=true)
                append!(supp[k][i], l*ones(UInt16, exp[j]))
            end
        end
    end
    if obj=="trace"
        supp[1], coe[1]=cyclic_canon(supp[1], coe[1])
    else
        supp[1], coe[1]=sym_canon(supp[1], coe[1])
    end
    lt[1]=length(supp[1])
    fbasis=get_ncbasis(n,d)
    gbasis=Vector{Vector{Vector{UInt16}}}(undef,m)
    tsupp=copy(supp[1])
    for i=1:m
        gbasis[i]=get_ncbasis(n,d-Int(ceil(maxdegree(pop[i+1])/2)))
        if obj=="trace"
            append!(tsupp, _cyclic_canon.(_sym_canon.(supp[i+1])))
        else
            append!(tsupp, _sym_canon.(supp[i+1]))
        end
    end
    if obj=="trace"
        append!(tsupp, [_cyclic_canon([fbasis[i][end:-1:1]; fbasis[i]]) for i=1:length(fbasis)])
    else
        append!(tsupp, [[fbasis[i][end:-1:1]; fbasis[i]] for i=1:length(fbasis)])
    end
    sort!(tsupp)
    unique!(tsupp)
    fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,_=get_nccblocks!(m,tsupp,supp[2:end],lt[2:end],fbasis,gbasis,TS=TS,obj=obj,QUIET=QUIET,merge=merge)
    if reducebasis==true&&obj=="eigen"
        gsupp=get_ncgsupp(m,lt,supp,gbasis,gblocks,gcl,gblocksize)
        psupp=copy(supp[1])
        push!(psupp, UInt16[])
        append!(psupp, gsupp)
        psupp=psupp[is_sym.(psupp)]
        fbasis,flag=reducebasis!(psupp,fbasis,fblocks,fcl,fblocksize)
        if flag==1
            tsupp=copy(supp[1])
            for i=1:m
                append!(tsupp, _sym_canon.(supp[i+1]))
            end
            append!(tsupp, [[fbasis[i][end:-1:1]; fbasis[i]] for i=1:length(fbasis)])
            sort!(tsupp)
            unique!(tsupp)
            fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,_=get_nccblocks!(m,tsupp,supp[2:end],lt[2:end],fbasis,gbasis,TS=TS,obj=obj,QUIET=QUIET,merge=merge)
        end
    end
    opt,tsupp=ncblockcpop(m,supp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,numeq=numeq,QUIET=QUIET,obj=obj)
    data=cdata_type(m,numeq,supp,coe,obj,lt,fbasis,gbasis,tsupp,ub,sizes,gblocks,gcl,gblocksize)
    return opt,data
end

function nctssos_first(supp::Vector{Vector{Vector{UInt16}}},coe::Vector{Vector{Float64}},n::Int64,d::Int64,dg::Vector{Int64};numeq=0,reducebasis=false,TS="block",obj="eigen",merge=false,QUIET=false)
    m=length(supp)-1
    if obj=="trace"
        supp[1], coe[1]=cyclic_canon(supp[1], coe[1])
    else
        supp[1], coe[1]=sym_canon(supp[1], coe[1])
    end
    lt=length.(supp)
    fbasis=get_ncbasis(n,d)
    gbasis=Vector{Vector{Vector{UInt16}}}(undef,m)
    tsupp=copy(supp[1])
    for i=1:m
        gbasis[i]=get_ncbasis(n,d-Int(ceil(dg[i]/2)))
        if obj=="trace"
            append!(tsupp, _cyclic_canon.(_sym_canon.(supp[i+1])))
        else
            append!(tsupp, _sym_canon.(supp[i+1]))
        end
    end
    if obj=="trace"
        append!(tsupp, [_cyclic_canon([fbasis[i][end:-1:1]; fbasis[i]]) for i=1:length(fbasis)])
    else
        append!(tsupp, [[fbasis[i][end:-1:1]; fbasis[i]] for i=1:length(fbasis)])
    end
    sort!(tsupp)
    unique!(tsupp)
    fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,_=get_nccblocks!(m,tsupp,supp[2:end],lt[2:end],fbasis,gbasis,TS=TS,obj=obj,QUIET=QUIET,merge=merge)
    if reducebasis==true&&obj=="eigen"
        gsupp=get_ncgsupp(m,lt,supp,gbasis,gblocks,gcl,gblocksize)
        psupp=copy(supp[1])
        push!(psupp, UInt16[])
        append!(psupp, gsupp)
        psupp=psupp[is_sym.(psupp)]
        fbasis,flag=reducebasis!(psupp,fbasis,fblocks,fcl,fblocksize)
        if flag==1
            tsupp=copy(supp[1])
            for i=1:m
                append!(tsupp, _sym_canon.(supp[i+1]))
            end
            append!(tsupp, [[fbasis[i][end:-1:1]; fbasis[i]] for i=1:length(fbasis)])
            sort!(tsupp)
            unique!(tsupp)
            fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,_=get_nccblocks!(m,tsupp,supp[2:end],lt[2:end],fbasis,gbasis,TS=TS,obj=obj,QUIET=QUIET,merge=merge)
        end
    end
    opt,tsupp=ncblockcpop(m,supp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,numeq=numeq,QUIET=QUIET,obj=obj)
    data=cdata_type(m,numeq,supp,coe,obj,lt,fbasis,gbasis,tsupp,ub,sizes,gblocks,gcl,gblocksize)
    return opt,data
end

function nctssos_higher!(data::cdata_type;TS="block",minimize=false,merge=false,QUIET=false)
    m=data.m
    numeq=data.numeq
    supp=data.supp
    coe=data.coe
    obj=data.obj
    lt=data.lt
    fbasis=data.fbasis
    gbasis=data.gbasis
    tsupp=data.tsupp
    ub=data.ub
    sizes=data.sizes
    gblocks=data.gblocks
    gcl=data.gcl
    gblocksize=data.gblocksize
    fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_nccblocks!(m,tsupp,supp[2:end],lt[2:end],fbasis,gbasis,gblocks=gblocks,gcl=gcl,gblocksize=gblocksize,ub=ub,sizes=sizes,TS=TS,obj=obj,QUIET=QUIET,merge=merge)
    opt=nothing
    if status==1
        opt,tsupp=ncblockcpop(m,supp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,numeq=numeq,QUIET=QUIET,obj=obj)
    end
    data.tsupp=tsupp
    data.ub=ub
    data.sizes=sizes
    data.gblocks=gblocks
    data.gcl=gcl
    data.gblocksize=gblocksize
    return opt,data
end

function get_ncgsupp(m,lt,supp,gbasis,gblocks,gcl,gblocksize)
    gsupp=Vector{UInt16}[]
    for k=1:m, i=1:gcl[k], j=1:gblocksize[k][i], r=j:gblocksize[k][i], s=1:lt[k+1]
        @inbounds bi=[gbasis[k][gblocks[k][i][j]][end:-1:1]; supp[k+1][s]; gbasis[k][gblocks[k][i][r]]]
        push!(gsupp, bi)
    end
    return gsupp
end

function reducebasis!(psupp,basis,blocks,cl,blocksize)
    init=0
    flag=0
    check=0
    while init==0||check>0
        init=1
        check=0
        for i=1:cl
            if blocksize[i]>1
                for j=1:blocksize[i], r=1:blocksize[i]
                    if j!=r
                        @inbounds bi = [basis[blocks[i][j]][end:-1:1]; basis[blocks[i][r]]]
                        if is_sym(bi)
                            push!(psupp, bi)
                        end
                    end
                end
            end
        end
        sort!(psupp)
        unique!(psupp)
        lpsupp=length(psupp)
        for i=1:cl
            lo=blocksize[i]
            indexb=[k for k=1:lo]
            j=1
            while lo>=j
                bi = [basis[blocks[i][indexb[j]]][end:-1:1]; basis[blocks[i][indexb[j]]]]
                Locb=ncbfind(psupp,lpsupp,bi)
                if Locb==0
                   check=1
                   flag=1
                   deleteat!(indexb, j)
                   lo=lo-1
                else
                   j=j+1
                end
            end
            blocks[i]=blocks[i][indexb]
            blocksize[i]=lo
        end
    end
    if flag==1
       indexb=blocks[1]
       for i=2:cl
           indexb=append!(indexb, blocks[i])
       end
       sort!(indexb)
       unique!(indexb)
       return basis[indexb],flag
    else
       return basis,flag
    end
end

function get_nccgraph(tsupp,supp,lt,basis;obj="eigen")
    lb=length(basis)
    ltsupp=length(tsupp)
    G=SimpleGraph(lb)
    for i = 1:lb, j = i+1:lb
        r=1
        while r<=lt
            bi=[basis[i][end:-1:1]; supp[r]; basis[j]]
            bi=_sym_canon(bi)
            if obj=="trace"
                bi=_cyclic_canon(bi)
            end
            if ncbfind(tsupp, ltsupp, bi)!=0
               break
            else
                r+=1
            end
        end
        if r<=lt
           add_edge!(G, i, j)
        end
    end
    return G
end

function get_nccblocks!(m,tsupp,gsupp,lt,fbasis,gbasis;gblocks=[],gcl=[],gblocksize=[],ub=[],sizes=[],TS="block",obj="eigen",minimize=false,QUIET=true,merge=false)
    if isempty(gblocks)
        gblocks=Vector{Vector{Vector{UInt16}}}(undef,m)
        gblocksize=Vector{Vector{UInt16}}(undef, m)
        gcl=Vector{UInt16}(undef,m)
    end
    if TS==false
        fblocksize=[length(fbasis)]
        fblocks=[[i for i=1:length(fbasis)]]
        fcl=1
        for k=1:m
            gblocks[k]=[[i for i=1:length(gbasis[k])]]
            gblocksize[k]=[length(gbasis[k])]
            gcl[k]=1
        end
        status=1
        nub=fblocksize
        nsizes=[1]
        if QUIET==false
            println("------------------------------------------------------")
            println("The sizes of blocks:\n$nub\n$nsizes")
            println("------------------------------------------------------")
        end
    elseif TS=="block"
        G=get_ncgraph(tsupp,fbasis,obj=obj)
        fblocks=connected_components(G)
        fblocksize=length.(fblocks)
        fcl=length(fblocksize)
        nub=unique(fblocksize)
        nsizes=[sum(fblocksize.== i) for i in nub]
        if isempty(ub)||nub!=ub||nsizes!=sizes
            status=1
            if QUIET==false
                println("------------------------------------------------------")
                println("The sizes of blocks:\n$nub\n$nsizes")
                println("------------------------------------------------------")
            end
            for k=1:m
                G=get_nccgraph(tsupp,gsupp[k],lt[k],gbasis[k],obj=obj)
                gblocks[k]=connected_components(G)
                gblocksize[k]=length.(gblocks[k])
                gcl[k]=length(gblocksize[k])
            end
        else
            status=0
            if QUIET==false
               println("No higher NCTSSOS hierarchy!")
            end
        end
    else
        G=get_ncgraph(tsupp,fbasis,obj=obj)
        fblocks,fcl,fblocksize=chordal_cliques!(G, method=TS, minimize=minimize)
        if merge==true
            fblocks,fcl,fblocksize=clique_merge!(fblocks,fcl,QUIET=true)
        end
        nub=unique(fblocksize)
        nsizes=[sum(fblocksize.== i) for i in nub]
        if isempty(ub)||nub!=ub||nsizes!=sizes
            status=1
            if QUIET==false
                println("------------------------------------------------------")
                println("The sizes of blocks:\n$nub\n$nsizes")
                println("------------------------------------------------------")
            end
            for k=1:m
                G=get_nccgraph(tsupp,gsupp[k],lt[k],gbasis[k],obj=obj)
                gblocks[k],gcl[k],gblocksize[k]=chordal_cliques!(G, method=TS, minimize=minimize)
                if merge==true
                    gblocks[k],gcl[k],gblocksize[k]=clique_merge!(gblocks[k],gcl[k],QUIET=true)
                end
            end
        else
            status=0
            if QUIET==false
               println("No higher NCTSSOS hierarchy!")
            end
        end
    end
    return fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,nub,nsizes,status
end

function ncblockcpop(m,supp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize;numeq=0,QUIET=true,obj="eigen")
    tsupp=Vector{UInt16}[]
    for i=1:fcl, j=1:fblocksize[i], r=j:fblocksize[i]
        @inbounds bi=[fbasis[fblocks[i][j]][end:-1:1]; fbasis[fblocks[i][r]]]
        @inbounds push!(tsupp, bi)
    end
    gsupp=get_ncgsupp(m,lt,supp,gbasis,gblocks,gcl,gblocksize)
    append!(tsupp, gsupp)
    tsupp=_sym_canon.(tsupp)
    if obj=="trace"
        tsupp=_cyclic_canon.(tsupp)
    end
    sort!(tsupp)
    unique!(tsupp)
    ltsupp=length(tsupp)
    model=Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons=[AffExpr(0) for i=1:ltsupp]
    pos=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, fcl)
    for i=1:fcl
        bs=fblocksize[i]
        if bs==1
           @inbounds pos[i]=@variable(model, lower_bound=0)
           @inbounds bi=[fbasis[fblocks[i][1]][end:-1:1]; fbasis[fblocks[i][1]]]
           if obj=="trace"
               bi=_cyclic_canon(bi)
           end
           Locb=ncbfind(tsupp,ltsupp,bi)
           @inbounds cons[Locb]+=pos[i]
        else
           @inbounds pos[i]=@variable(model, [1:bs, 1:bs], PSD)
           for j=1:bs, r=j:bs
               @inbounds bi=[fbasis[fblocks[i][j]][end:-1:1]; fbasis[fblocks[i][r]]]
               bi=_sym_canon(bi)
               if obj=="trace"
                   bi=_cyclic_canon(bi)
               end
               Locb=ncbfind(tsupp,ltsupp,bi)
               if j==r
                   @inbounds cons[Locb]+=pos[i][j,r]
               else
                   @inbounds cons[Locb]+=2*pos[i][j,r]
               end
           end
        end
    end
    gpos=Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, m)
    for k=1:m
        gpos[k]=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, gcl[k])
        for i=1:gcl[k]
            bs=gblocksize[k][i]
            if bs==1
                if k<=m-numeq
                    gpos[k][i]=@variable(model, lower_bound=0)
                else
                    gpos[k][i]=@variable(model)
                end
                for s=1:lt[k+1]
                    @inbounds bi=[gbasis[k][gblocks[k][i][1]][end:-1:1]; supp[k+1][s]; gbasis[k][gblocks[k][i][1]]]
                    bi=_sym_canon(bi)
                    if obj=="trace"
                        bi=_cyclic_canon(bi)
                    end
                    Locb=ncbfind(tsupp,ltsupp,bi)
                    @inbounds cons[Locb]+=coe[k+1][s]*gpos[k][i]
                end
            else
                if k<=m-numeq
                   gpos[k][i]=@variable(model, [1:bs, 1:bs], PSD)
                else
                   gpos[k][i]=@variable(model, [1:bs, 1:bs], Symmetric)
                end
                for j=1:bs, r=j:bs, s=1:lt[k+1]
                    @inbounds bi=[gbasis[k][gblocks[k][i][j]][end:-1:1]; supp[k+1][s]; gbasis[k][gblocks[k][i][r]]]
                    bi=_sym_canon(bi)
                    if obj=="trace"
                        bi=_cyclic_canon(bi)
                    end
                    Locb=ncbfind(tsupp,ltsupp,bi)
                    if j==r
                        @inbounds cons[Locb]+=coe[k+1][s]*gpos[k][i][j,r]
                    else
                        @inbounds cons[Locb]+=2*coe[k+1][s]*gpos[k][i][j,r]
                    end
                end
            end
        end
    end
    bc=zeros(ltsupp)
    for i=1:lt[1]
        Locb=ncbfind(tsupp,ltsupp,supp[1][i])
        if Locb==0
           @error "The monomial basis is not enough!"
           return nothing,nothing
        else
           bc[Locb]=coe[1][i]
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
