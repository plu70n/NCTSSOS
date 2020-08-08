mutable struct cdata_type
    m
    ssupp
    coe
    obj
    lt
    fbasis
    gbasis
    fsupp
    ub
    sizes
    numeq
    gblocks
    gcl
    gblocksize
end

function ncblockcpop_first(pop,x,d;numeq=0,reducebasis=false,TS="block",obj="eigen",merge=false,QUIET=false)
    n=length(x)
    if obj=="trace"
        pop[1]=cyclic_canon(pop[1], x)
    end
    m=length(pop)-1
    dg=zeros(Int,m)
    coe=Vector{Vector{Float64}}(undef, m+1)
    ssupp=Vector{Vector{Vector{UInt16}}}(undef, m+1)
    lt=Vector{UInt16}(undef, m+1)
    for k=1:m+1
        mon=monomials(pop[k])
        coe[k]=coefficients(pop[k])
        lt[k]=length(mon)
        ssupp[k]=[UInt16[] for i=1:lt[k]]
        for i=1:lt[k]
            vars=variables(mon[i])
            exp=exponents(mon[i])
            ind=[exp[j]!=0 for j=1:length(exp)]
            vars=vars[ind]
            exp=exp[ind]
            for j=1:length(vars)
                l=ncbfind(x, n, vars[j], rev=true)
                append!(ssupp[k][i], l*ones(UInt16, exp[j]))
            end
        end
    end
    supp=ssupp[1]
    for i=2:m+1
        dg[i-1]=maxdegree(pop[i])
        append!(supp, ssupp[i])
    end
    unique!(supp)
    fbasis=get_ncbasis(n,d)
    gbasis=Vector{Vector{Vector{UInt16}}}(undef,m)
    for k=1:m
        gbasis[k]=get_ncbasis(n,d-Int(ceil(dg[k]/2)))
    end
    fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,_=get_nccblocks!(m,supp,ssupp,lt,fbasis,gbasis,TS=TS,obj=obj,QUIET=QUIET,merge=merge)
    if reducebasis==true
        gsupp=get_ncgsupp(m,lt,ssupp,gbasis,gblocks,gcl,gblocksize)
        tsupp=copy(ssupp[1])
        push!(tsupp, UInt16[])
        append!(tsupp, gsupp)
        fbasis,flag=reducebasis!(tsupp,fbasis,fblocks,fcl,fblocksize)
        if flag==1
            fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,_=get_nccblocks!(m,supp,ssupp,lt,fbasis,gbasis,TS=TS,obj=obj,QUIET=QUIET,merge=merge)
        end
    end
    opt,fsupp=ncblockcpop(m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,numeq=numeq,QUIET=QUIET,obj=obj)
    data=cdata_type(m,ssupp,coe,obj,lt,fbasis,gbasis,fsupp,ub,sizes,numeq,gblocks,gcl,gblocksize)
    return opt,data
end

function ncblockcpop_higher!(data;TS="block",minimize=false,merge=false,QUIET=false)
    m=data.m
    ssupp=data.ssupp
    coe=data.coe
    obj=data.obj
    lt=data.lt
    fbasis=data.fbasis
    gbasis=data.gbasis
    fsupp=data.fsupp
    ub=data.ub
    sizes=data.sizes
    numeq=data.numeq
    gblocks=data.gblocks
    gcl=data.gcl
    gblocksize=data.gblocksize
    fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,ub,sizes,status=get_nccblocks!(m,fsupp,ssupp,lt,fbasis,gbasis,gblocks=gblocks,gcl=gcl,gblocksize=gblocksize,ub=ub,sizes=sizes,TS=TS,obj=obj,QUIET=QUIET,merge=merge)
    opt=nothing
    if status==1
        opt,fsupp=ncblockcpop(m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,numeq=numeq,QUIET=QUIET,obj=obj)
    end
    data.fsupp=fsupp
    data.ub=ub
    data.sizes=sizes
    data.gblocks=gblocks
    data.gcl=gcl
    data.gblocksize=gblocksize
    return opt,data
end

function get_ncgsupp(m,lt,ssupp,gbasis,gblocks,gcl,gblocksize)
    gsupp=Vector{Vector{UInt16}}(undef, sum(lt[k+1]*sum(gblocksize[k].^2) for k=1:m))
    l=1
    for k=1:m, i=1:gcl[k], j=1:gblocksize[k][i], r=1:gblocksize[k][i], s=1:lt[k+1]
        @inbounds gsupp[l]=[gbasis[k][gblocks[k][i][j]][end:-1:1]; ssupp[k+1][s]; gbasis[k][gblocks[k][i][r]]]
        l+=1
    end
    return gsupp
end

function reducebasis!(tsupp,basis,blocks,cl,blocksize)
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
                        tsupp=push!(tsupp, bi)
                    end
                end
            end
        end
        sort!(tsupp)
        unique!(tsupp)
        ltsupp=length(tsupp)
        for i=1:cl
            lo=blocksize[i]
            indexb=[k for k=1:lo]
            j=1
            while lo>=j
                bi = [basis[blocks[i][indexb[j]]][end:-1:1]; basis[blocks[i][indexb[j]]]]
                Locb=ncbfind(tsupp,ltsupp,bi)
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

function get_nccgraph(tsupp,ssupp,lt,basis;obj="eigen")
    lb=length(basis)
    ltsupp=length(tsupp)
    G=SimpleGraph(lb)
    for i = 1:lb, j = i+1:lb
        r=1
        while r<=lt
            bi=[basis[i][end:-1:1]; ssupp[r]; basis[j]]
            if obj=="trace"
                bi=_cyclic_canon(bi)
            end
            if ncbfind(tsupp,ltsupp,bi)!=0
               break
            else
                r+=1
            end
        end
        if r<=lt
           add_edge!(G,i,j)
        end
    end
    return G
end

function get_nccblocks!(m,supp,ssupp,lt,fbasis,gbasis;gblocks=[],gcl=[],gblocksize=[],ub=[],sizes=[],TS="block",obj="eigen",minimize=false,QUIET=true,merge=false)
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
        G,tsupp=get_ncgraph(supp,fbasis,obj=obj)
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
                G=get_nccgraph(tsupp,ssupp[k+1],lt[k+1],gbasis[k],obj=obj)
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
        G,tsupp=get_ncgraph(supp,fbasis,obj=obj)
        fblocks,fcl,fblocksize=chordal_cliques!(G, method=TS, minimize=minimize)
        if merge==true
            fblocks,fcl,fblocksize=clique_merge!(fblocks,fcl,QUIET=true)
        end
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
                G=get_nccgraph(tsupp,ssupp[k+1],lt[k+1],gbasis[k],obj=obj)
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

function ncblockcpop(m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize;numeq=0,QUIET=true,obj="eigen")
    fsupp=Vector{Vector{UInt16}}(undef, sum(fblocksize.^2))
    k=1
    for i=1:fcl, j=1:fblocksize[i], r=1:fblocksize[i]
        @inbounds bi=[fbasis[fblocks[i][j]][end:-1:1]; fbasis[fblocks[i][r]]]
        @inbounds fsupp[k]=bi
        k+=1
    end
    gsupp=get_ncgsupp(m,lt,ssupp,gbasis,gblocks,gcl,gblocksize)
    supp1=append!(fsupp, gsupp)
    if obj=="trace"
        supp1=_cyclic_canon.(supp1)
    end
    supp1=sort!(supp1)
    supp1=unique!(supp1)
    lsupp1=length(supp1)
    model=Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons=[AffExpr(0) for i=1:lsupp1]
    pos=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, fcl)
    for i=1:fcl
        bs=fblocksize[i]
        if bs==1
           @inbounds pos[i]=@variable(model, lower_bound=0)
           @inbounds bi=[fbasis[fblocks[i][1]][end:-1:1]; fbasis[fblocks[i][1]]]
           if obj=="trace"
               bi=_cyclic_canon(bi)
           end
           Locb=ncbfind(supp1,lsupp1,bi)
           @inbounds cons[Locb]+=pos[i]
        else
           @inbounds pos[i]=@variable(model, [1:bs, 1:bs], PSD)
           for j=1:bs, r=1:bs
               @inbounds bi=[fbasis[fblocks[i][j]][end:-1:1]; fbasis[fblocks[i][r]]]
               if obj=="trace"
                   bi=_cyclic_canon(bi)
               end
               Locb=ncbfind(supp1,lsupp1,bi)
               @inbounds cons[Locb]+=pos[i][j,r]
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
                    @inbounds bi=[gbasis[k][gblocks[k][i][1]][end:-1:1]; ssupp[k+1][s]; gbasis[k][gblocks[k][i][1]]]
                    if obj=="trace"
                        bi=_cyclic_canon(bi)
                    end
                    Locb=ncbfind(supp1,lsupp1,bi)
                    @inbounds cons[Locb]+=coe[k+1][s]*gpos[k][i]
                end
            else
                if k<=m-numeq
                   gpos[k][i]=@variable(model, [1:bs, 1:bs], PSD)
                else
                   gpos[k][i]=@variable(model, [1:bs, 1:bs], Symmetric)
                end
                for j=1:bs, r=1:bs, s=1:lt[k+1]
                    @inbounds bi=[gbasis[k][gblocks[k][i][j]][end:-1:1]; ssupp[k+1][s]; gbasis[k][gblocks[k][i][r]]]
                    if obj=="trace"
                        bi=_cyclic_canon(bi)
                    end
                    Locb=ncbfind(supp1,lsupp1,bi)
                    @inbounds cons[Locb]+=coe[k+1][s]*gpos[k][i][j,r]
                end
            end
        end
    end
    bc=zeros(lsupp1)
    for i=1:lt[1]
        Locb=ncbfind(supp1,lsupp1,ssupp[1][i])
        if Locb==0
           @error "The monomial basis is not enough!"
           return nothing,nothing
        else
           bc[Locb]=coe[1][i]
       end
    end
    @variable(model, lower)
    cons[1]+=lower
    @constraint(model, con[i=1:lsupp1], cons[i]==bc[i])
    @objective(model, Max, lower)
    optimize!(model)
    status=termination_status(model)
    if  status==MOI.OPTIMAL
        objv=objective_value(model)
        println("optimum = $objv")
    else
        objv=objective_value(model)
        println("termination status: $status")
        sstatus=primal_status(model)
        println("solution status: $sstatus")
        println("optimum = $objv")
    end
    return objv,fsupp
end
