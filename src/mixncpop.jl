mutable struct mdata_type
    n::Int
    m::Int
    d::Int
    dg::Vector{Int}
    supp::Vector{Vector{Vector{UInt16}}}
    coe::Vector{Vector{Float64}}
    obj
    numeq::Int
    tsupp::Vector{Vector{UInt16}}
    lt::Vector{Vector{UInt16}}
    fbasis::Vector{Vector{Vector{UInt16}}}
    gbasis::Vector{Vector{Vector{Vector{UInt16}}}}
    cql::Int
    cliques::Vector{Vector{UInt16}}
    cliquesize::Vector{Int}
    J::Vector{Vector{UInt16}}
    ncc::Vector{UInt16}
    blocks::Vector{Vector{Vector{Vector{UInt16}}}}
    cl::Vector{Vector{UInt16}}
    blocksize::Vector{Vector{Vector{Int}}}
    ub::Vector{Vector{UInt16}}
    sizes::Vector{Vector{UInt16}}
end

mutable struct mudata_type
    n::Int
    d::Int
    supp::Vector{Vector{UInt16}}
    coe::Vector{Float64}
    obj
    tsupp::Vector{Vector{UInt16}}
    basis::Vector{Vector{Vector{UInt16}}}
    cql::Int
    cliques::Vector{Vector{UInt16}}
    cliquesize::Vector{Int}
    blocks::Vector{Vector{Vector{UInt16}}}
    cl::Vector{UInt16}
    blocksize::Vector{Vector{Int}}
    ub::Vector{Vector{UInt16}}
    sizes::Vector{Vector{UInt16}}
end

function cs_nctssos_first(f,x;d=0,CS="MD",minimize=false,TS="block",QUIET=false,obj="eigen",solve=true)
    n=length(x)
    mon=monomials(f)
    coe=coefficients(f)
    supp=[UInt16[] for i=1:length(mon)]
    for i=1:length(mon)
        ind=mon[i].z .>0
        vars=mon[i].vars[ind]
        exp=mon[i].z[ind]
        for j=1:length(vars)
            l=ncbfind(x, n, vars[j], rev=true)
            append!(supp[i], l*ones(UInt16, exp[j]))
        end
    end
    if obj=="trace"
        supp, coe=cyclic_canon(supp, coe)
    else
        supp, coe=sym_canon(supp, coe)
    end
    if d==0
        d=ceil(Int, maxdegree(f)/2)
    end
    cliques,cql,cliquesize=clique_decomp(n,supp)
    blocks,cl,blocksize,ub,sizes,basis,_=get_blocks_mix(d,supp,cliques,cql,cliquesize,TS=TS,obj=obj)
    opt,tsupp=blockupop_mix(n,d,supp,coe,basis,cliques,cql,cliquesize,blocks,cl,blocksize,obj=obj,solve=solve)
    data=mudata_type(n,d,supp,coe,obj,tsupp,basis,cql,cliques,cliquesize,blocks,cl,blocksize,ub,sizes)
    return opt,data
end

function cs_nctssos_first(supp::Vector{Vector{UInt16}},coe::Vector{Float64},n::Int;d=0,CS="MD",minimize=false,TS="block",QUIET=false,obj="eigen",solve=true)
    if obj=="trace"
        supp, coe=cyclic_canon(supp, coe)
    else
        supp, coe=sym_canon(supp, coe)
    end
    cliques,cql,cliquesize=clique_decomp(n,supp)
    blocks,cl,blocksize,ub,sizes,basis,_=get_blocks_mix(d,supp,cliques,cql,cliquesize,TS=TS,obj=obj)
    opt,tsupp=blockupop_mix(n,d,supp,coe,basis,cliques,cql,cliquesize,blocks,cl,blocksize,obj=obj,solve=solve)
    data=mudata_type(n,d,supp,coe,obj,tsupp,basis,cql,cliques,cliquesize,blocks,cl,blocksize,ub,sizes)
    return opt,data
end

function cs_nctssos_higher!(data::mudata_type;TS="block",merge=false,QUIET=false,solve=true)
    n=data.n
    d=data.d
    supp=data.supp
    coe=data.coe
    obj=data.obj
    tsupp=data.tsupp
    basis=data.basis
    cql=data.cql
    cliques=data.cliques
    cliquesize=data.cliquesize
    blocks=data.blocks
    cl=data.cl
    blocksize=data.blocksize
    ub=data.ub
    sizes=data.sizes
    blocks,cl,blocksize,ub,sizes,basis,status=get_blocks_mix(d,supp,cliques,cql,cliquesize,basis=basis,ub=ub,sizes=sizes,TS=TS,merge=merge,obj=obj)
    if status==1
        opt,tsupp=blockupop_mix(n,d,supp,coe,basis,cliques,cql,cliquesize,blocks,cl,blocksize,obj=obj,solve=solve)
    else
        opt=nothing
        println("No higher CS-NCTSSOS hierarchy!")
    end
    data.tsupp=tsupp
    data.blocks=blocks
    data.cl=cl
    data.blocksize=blocksize
    data.ub=ub
    data.sizes=sizes
    return opt,data
end

function cs_nctssos_first(pop,x,d;numeq=0,CS="MD",minimize=false,assign="min",TS="block",QUIET=false,obj="eigen",solve=true)
    n=length(x)
    m=length(pop)-1
    coe=Vector{Vector{Float64}}(undef, m+1)
    supp=Vector{Vector{Vector{UInt16}}}(undef, m+1)
    for k=1:m+1
        mon=monomials(pop[k])
        coe[k]=coefficients(pop[k])
        supp[k]=[UInt16[] for i=1:length(mon)]
        for i=1:length(mon)
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
    dg=zeros(Int,m)
    for i=1:m
        dg[i]=maxdegree(pop[i+1])
    end
    cliques,cql,cliquesize=clique_decomp(n,m,d,dg,supp,alg=CS,minimize=minimize)
    J,ncc=assign_constraint(m,supp,cliques,cql,cliquesize,assign=assign)
    blocks,cl,blocksize,ub,sizes,lt,fbasis,gbasis,status=get_cblocks_mix!(d,dg,J,m,supp,cliques,cql,cliquesize,TS=TS,obj=obj)
    opt,tsupp=blockcpop_mix(n,m,d,dg,supp,coe,fbasis,gbasis,cliques,cql,cliquesize,J,ncc,blocks,cl,blocksize,numeq=numeq,QUIET=QUIET,obj=obj,solve=solve)
    data=mdata_type(n,m,d,dg,supp,coe,obj,numeq,tsupp,lt,fbasis,gbasis,cql,cliques,cliquesize,J,ncc,blocks,cl,blocksize,ub,sizes)
    return opt,data
end

function cs_nctssos_first(supp::Vector{Vector{Vector{UInt16}}},coe::Vector{Vector{Float64}},n::Int,d::Int,dg::Vector{Int};numeq=0,CS="MD",minimize=false,assign="min",TS="block",QUIET=false,obj="eigen",solve=true)
    m=length(supp)-1
    if obj=="trace"
        supp[1], coe[1]=cyclic_canon(supp[1], coe[1])
    else
        supp[1], coe[1]=sym_canon(supp[1], coe[1])
    end
    cliques,cql,cliquesize=clique_decomp(n,m,d,dg,supp,alg=CS,minimize=minimize)
    J,ncc=assign_constraint(m,supp,cliques,cql,cliquesize,assign=assign)
    blocks,cl,blocksize,ub,sizes,lt,fbasis,gbasis,status=get_cblocks_mix!(d,dg,J,m,supp,cliques,cql,cliquesize,TS=TS,obj=obj)
    opt,tsupp=blockcpop_mix(n,m,d,dg,supp,coe,fbasis,gbasis,cliques,cql,cliquesize,J,ncc,blocks,cl,blocksize,numeq=numeq,QUIET=QUIET,obj=obj,solve=solve)
    data=mdata_type(n,m,d,dg,supp,coe,obj,numeq,tsupp,lt,fbasis,gbasis,cql,cliques,cliquesize,J,ncc,blocks,cl,blocksize,ub,sizes)
    return opt,data
end

function cs_nctssos_higher!(data::mdata_type;TS="block",QUIET=false,solve=true)
    n=data.n
    m=data.m
    d=data.d
    dg=data.dg
    supp=data.supp
    coe=data.coe
    obj=data.obj
    numeq=data.numeq
    tsupp=data.tsupp
    lt=data.lt
    fbasis=data.fbasis
    gbasis=data.gbasis
    cql=data.cql
    cliques=data.cliques
    cliquesize=data.cliquesize
    J=data.J
    ncc=data.ncc
    blocks=data.blocks
    cl=data.cl
    blocksize=data.blocksize
    ub=data.ub
    sizes=data.sizes
    blocks,cl,blocksize,ub,sizes,lt,fbasis,gbasis,status=get_cblocks_mix!(d,dg,J,m,supp,cliques,cql,cliquesize,tsupp=tsupp,lt=lt,fbasis=fbasis,gbasis=gbasis,blocks=blocks,cl=cl,blocksize=blocksize,ub=ub,sizes=sizes,TS=TS,obj=obj)
    if status==1
        opt,tsupp=blockcpop_mix(n,m,d,dg,supp,coe,fbasis,gbasis,cliques,cql,cliquesize,J,ncc,blocks,cl,blocksize,numeq=numeq,QUIET=QUIET,obj=obj,solve=solve)
    else
        opt=nothing
        println("No higher CS-NCTSSOS hierarchy!")
    end
    data.tsupp=tsupp
    data.blocks=blocks
    data.cl=cl
    data.blocksize=blocksize
    data.ub=ub
    data.sizes=sizes
    return opt,data
end

function blockupop_mix(n,d,supp,coe,basis,cliques,cql,cliquesize,blocks,cl,blocksize;QUIET=false,obj="eigen",solve=true)
    tsupp=Vector{UInt16}[]
    for i=1:cql, j=1:cl[i], k=1:blocksize[i][j], r=k:blocksize[i][j]
        @inbounds bi=[basis[i][blocks[i][j][k]][end:-1:1]; basis[i][blocks[i][j][r]]]
        push!(tsupp, bi)
    end
    tsupp=_sym_canon.(tsupp)
    if obj=="trace"
        tsupp=_cyclic_canon.(tsupp)
    end
    # if nx>0
    #     tsupp=comm.(tsupp, nx)
    #     proj!.(tsupp)
    # end
    sort!(tsupp)
    unique!(tsupp)
    objv=nothing
    if solve==true
        ltsupp=length(tsupp)
        model=Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        cons=[AffExpr(0) for i=1:ltsupp]
        pos1=Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, cql)
        for i=1:cql
            pos1[i]=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, cl[i])
            for k=1:cl[i]
                if blocksize[i][k]==1
                   pos1[i][k]=@variable(model, lower_bound=0)
                   @inbounds bi=[basis[i][blocks[i][k][1]][end:-1:1]; basis[i][blocks[i][k][1]]]
                   if obj=="trace"
                       bi=_cyclic_canon(bi)
                   end
                   # if nx>0
                   #     bi=comm(bi, nx)
                   #     proj!(bi)
                   # end
                   Locb=ncbfind(tsupp,ltsupp,bi)
                   @inbounds cons[Locb]+=pos1[i][k]
                else
                   pos1[i][k]=@variable(model, [1:blocksize[i][k], 1:blocksize[i][k]], PSD)
                   for j=1:blocksize[i][k], r=j:blocksize[i][k]
                       @inbounds ind1=blocks[i][k][j]
                       @inbounds ind2=blocks[i][k][r]
                       @inbounds bi=[basis[i][ind1][end:-1:1]; basis[i][ind2]]
                       bi=_sym_canon(bi)
                       if obj=="trace"
                           bi=_cyclic_canon(bi)
                       end
                       # if nx>0
                       #     bi=comm(bi, nx)
                       #     proj!(bi)
                       # end
                       Locb=ncbfind(tsupp,ltsupp,bi)
                       if j==r
                           @inbounds cons[Locb]+=pos1[i][k][j,r]
                       else
                           @inbounds cons[Locb]+=2*pos1[i][k][j,r]
                       end
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
    end
    return objv,tsupp
end

function blockcpop_mix(n,m,d,dg,supp,coe,fbasis,gbasis,cliques,cql,cliquesize,J,ncc,blocks,cl,blocksize;numeq=0,QUIET=false,obj="eigen",solve=true)
    tsupp=Vector{UInt16}[]
    for i=1:cql
        for j=1:cl[i][1], k=1:blocksize[i][1][j], r=k:blocksize[i][1][j]
            @inbounds bi=[fbasis[i][blocks[i][1][j][k]][end:-1:1]; fbasis[i][blocks[i][1][j][r]]]
            push!(tsupp, bi)
        end
        for (j, w) in enumerate(J[i])
            for l=1:cl[i][j+1], t=1:blocksize[i][j+1][l], r=t:blocksize[i][j+1][l], s=1:length(supp[w+1])
                ind1=blocks[i][j+1][l][t]
                ind2=blocks[i][j+1][l][r]
                @inbounds bi=[gbasis[i][j][ind1][end:-1:1]; supp[w+1][s]; gbasis[i][j][ind2]]
                push!(tsupp, bi)
            end
        end
    end
    for i ∈ ncc
        append!(tsupp, supp[i+1])
    end
    tsupp=_sym_canon.(tsupp)
    if obj=="trace"
        tsupp=_cyclic_canon.(tsupp)
    end
    sort!(tsupp)
    unique!(tsupp)
    objv=nothing
    if solve==true
        ltsupp=length(tsupp)
        model=Model(optimizer_with_attributes(Mosek.Optimizer))
        set_optimizer_attribute(model, MOI.Silent(), QUIET)
        cons=[AffExpr(0) for i=1:ltsupp]
        for i=1:cql, l=1:cl[i][1]
            if blocksize[i][1][l]==1
               @inbounds pos=@variable(model, lower_bound=0)
               @inbounds bi=[fbasis[i][blocks[i][1][l][1]][end:-1:1]; fbasis[i][blocks[i][1][l][1]]]
               if obj=="trace"
                   bi=_cyclic_canon(bi)
               end
               Locb=ncbfind(tsupp,ltsupp,bi)
               @inbounds cons[Locb]+=pos
            else
               @inbounds bs=blocksize[i][1][l]
               @inbounds pos=@variable(model, [1:bs, 1:bs], PSD)
               for t=1:bs, r=t:bs
                   @inbounds ind1=blocks[i][1][l][t]
                   @inbounds ind2=blocks[i][1][l][r]
                   @inbounds bi=[fbasis[i][ind1][end:-1:1]; fbasis[i][ind2]]
                   bi=_sym_canon(bi)
                   if obj=="trace"
                       bi=_cyclic_canon(bi)
                   end
                   Locb=ncbfind(tsupp, ltsupp, bi)
                   if t==r
                      @inbounds cons[Locb]+=pos[t,r]
                   else
                      @inbounds cons[Locb]+=2*pos[t,r]
                   end
               end
            end
        end
        for k=1:length(ncc)
            i=ncc[k]
            if i<=m-numeq
                pos=@variable(model, lower_bound=0)
            else
                pos=@variable(model)
            end
            for j=1:length(supp[i+1])
                bi=_sym_canon(supp[i+1][j])
                if obj=="trace"
                    bi=_cyclic_canon(bi)
                end
                Locb=ncbfind(tsupp,ltsupp,bi)
                cons[Locb]+=coe[i+1][j]*pos
            end
        end
        for i=1:cql, (j, w) in enumerate(J[i])
            for l=1:cl[i][j+1]
                bs=blocksize[i][j+1][l]
                if bs==1
                    if j<=m-numeq
                        pos=@variable(model, lower_bound=0)
                    else
                        pos=@variable(model)
                    end
                    ind1=blocks[i][j+1][l][1]
                    for s=1:length(supp[w+1])
                        @inbounds bi=[gbasis[i][j][ind1][end:-1:1]; supp[w+1][s]; gbasis[i][j][ind1]]
                        bi=_sym_canon(bi)
                        if obj=="trace"
                            bi=_cyclic_canon(bi)
                        end
                        Locb=ncbfind(tsupp,ltsupp,bi)
                        @inbounds cons[Locb]+=coe[w+1][s]*pos
                    end
                else
                    if j<=m-numeq
                        pos=@variable(model, [1:bs, 1:bs], PSD)
                    else
                        pos=@variable(model, [1:bs, 1:bs], Symmetric)
                    end
                    for t=1:bs, r=t:bs
                        ind1=blocks[i][j+1][l][t]
                        ind2=blocks[i][j+1][l][r]
                        for s=1:length(supp[w+1])
                            @inbounds bi=[gbasis[i][j][ind1][end:-1:1]; supp[w+1][s]; gbasis[i][j][ind2]]
                            bi=_sym_canon(bi)
                            if obj=="trace"
                                bi=_cyclic_canon(bi)
                            end
                            Locb=ncbfind(tsupp,ltsupp,bi)
                            if t==r
                                @inbounds cons[Locb]+=coe[w+1][s]*pos[t,r]
                            else
                                @inbounds cons[Locb]+=2*coe[w+1][s]*pos[t,r]
                            end
                        end
                    end
                end
            end
        end
        bc=zeros(ltsupp)
        for i=1:length(supp[1])
            Locb=ncbfind(tsupp,ltsupp,supp[1][i])
            if Locb==0
               @error "The monomial basis is not enough!"
               return nothing,nothing
            else
               bc[Locb]=coe[1][i]
            end
        end
        @variable(model, lower)
        @constraint(model, cons[2:end].==bc[2:end])
        @constraint(model, cons[1]+lower==bc[1])
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
    end
    return objv,tsupp
end

function get_blocks_mix(d,supp,cliques,cql,cliquesize;basis=[],ub=[],sizes=[],TS="block",merge=false,obj="eigen")
    blocks=Vector{Vector{Vector{UInt16}}}(undef,cql)
    cl=Vector{UInt16}(undef,cql)
    blocksize=Vector{Vector{Int}}(undef,cql)
    status=ones(UInt8,cql)
    if isempty(basis)
        ub=Vector{Vector{UInt16}}(undef,cql)
        sizes=Vector{Vector{UInt16}}(undef,cql)
        basis=Vector{Vector{Vector{UInt16}}}(undef,cql)
        flag=1
    else
        flag=0
    end
    for i=1:cql
        nvar=cliquesize[i]
        ind=[issubset(supp[j], cliques[i]) for j=1:length(supp)]
        tsupp=copy(supp[ind])
        if flag==1
            basis[i]=get_ncbasis(nvar, d, ind=cliques[i])
            # if nx>0
            #      basis[i]=basis[i][is_basis.(basis[i], nx)]
            # end
            if obj=="trace"
                append!(tsupp, [_cyclic_canon([basis[i][k][end:-1:1]; basis[i][k]]) for k=1:length(basis[i])])
            else
                append!(tsupp, [[basis[i][k][end:-1:1]; basis[i][k]] for k=1:length(basis[i])])
            end
            # if nx>0
            #     tsupp=comm.(tsupp, nx)
            #     proj!.(tsupp)
            # end
            sort!(tsupp)
            unique!(tsupp)
            blocks[i],cl[i],blocksize[i],ub[i],sizes[i],status[i]=get_ncblocks(tsupp,basis[i],TS=TS,obj=obj,QUIET=true,merge=merge)
        else
            blocks[i],cl[i],blocksize[i],ub[i],sizes[i],status[i]=get_ncblocks(tsupp,basis[i],ub=ub[i],sizes=sizes[i],TS=TS,obj=obj,QUIET=true,merge=merge)
        end
    end
    return blocks,cl,blocksize,ub,sizes,basis,maximum(status)
end

function get_cblocks_mix!(d,dg,J,m,supp,cliques,cql,cliquesize;tsupp=[],lt=[],fbasis=[],gbasis=[],blocks=[],cl=[],blocksize=[],ub=[],sizes=[],TS="block",merge=false,obj="eigen")
    if isempty(fbasis)
        blocks=Vector{Vector{Vector{Vector{UInt16}}}}(undef,cql)
        cl=Vector{Vector{UInt16}}(undef,cql)
        blocksize=Vector{Vector{Vector{Int}}}(undef,cql)
        ub=Vector{Vector{UInt16}}(undef,cql)
        sizes=Vector{Vector{UInt16}}(undef,cql)
        lt=Vector{Vector{UInt16}}(undef,cql)
        fbasis=Vector{Vector{Vector{UInt16}}}(undef,cql)
        gbasis=Vector{Vector{Vector{Vector{UInt16}}}}(undef,cql)
        tsupp=copy(supp[1])
        for i=2:m+1
            append!(tsupp,  _sym_canon.(supp[i]))
        end
        sort!(tsupp)
        unique!(tsupp)
        flag=1
    else
        flag=0
    end
    status=ones(UInt8, cql)
    for i=1:cql
        lc=length(J[i])
        nvar=cliquesize[i]
        ind=[issubset(tsupp[j], cliques[i]) for j=1:length(tsupp)]
        fsupp=copy(tsupp[ind])
        if flag==1
            fbasis[i]=get_ncbasis(cliquesize[i], d, ind=cliques[i])
            gbasis[i]=Vector{Vector{Vector{UInt16}}}(undef, lc)
            lt[i]=length.(supp[J[i].+1])
            for s=1:lc
                gbasis[i][s]=get_ncbasis(nvar, d-ceil(Int, dg[J[i][s]]/2), ind=cliques[i])
            end
            blocks[i]=Vector{Vector{Vector{UInt16}}}(undef, lc+1)
            cl[i]=Vector{UInt16}(undef, lc+1)
            blocksize[i]=Vector{Vector{Int}}(undef, lc+1)
            ub[i]=Vector{UInt16}(undef, lc+1)
            sizes[i]=Vector{UInt16}(undef, lc+1)
            if obj=="trace"
                append!(fsupp, [_cyclic_canon([fbasis[i][k][end:-1:1]; fbasis[i][k]]) for k=1:length(fbasis[i])])
            else
                append!(fsupp, [[fbasis[i][k][end:-1:1]; fbasis[i][k]] for k=1:length(fbasis[i])])
            end
            sort!(fsupp)
            unique!(fsupp)
            blocks[i][1],cl[i][1],blocksize[i][1],blocks[i][2:end],cl[i][2:end],blocksize[i][2:end],ub[i],sizes[i],status[i]=get_nccblocks!(lc,fsupp,supp[J[i].+1],lt[i],fbasis[i],gbasis[i],TS=TS,QUIET=true,merge=merge,obj=obj)
        else
            blocks[i][1],cl[i][1],blocksize[i][1],blocks[i][2:end],cl[i][2:end],blocksize[i][2:end],ub[i],sizes[i],status[i]=get_nccblocks!(lc,fsupp,supp[J[i].+1],lt[i],fbasis[i],gbasis[i],gblocks=blocks[i][2:end],gcl=cl[i][2:end],gblocksize=blocksize[i][2:end],ub=ub[i],sizes=sizes[i],TS=TS,QUIET=true,merge=merge,obj=obj)
        end
    end
    return blocks,cl,blocksize,ub,sizes,lt,fbasis,gbasis,maximum(status)
end

function assign_constraint(m,supp,cliques,cql,cliquesize;assign="first")
    J=[UInt16[] for i=1:cql]
    ncc=UInt16[]
    for i=2:m+1
        rind=copy(supp[i][1])
        for j=2:length(supp[i])
            append!(rind, supp[i][j])
        end
        rind=unique(rind)
        if assign=="first"
            ind=findfirst(k->issubset(rind, cliques[k]), 1:cql)
            if ind!=nothing
                push!(J[ind], i-1)
            else
                push!(ncc, i-1)
            end
        else
            temp=UInt16[]
            for j=1:cql
                if issubset(rind, cliques[j])
                    push!(temp,j)
                end
            end
            if !isempty(temp)
                if assign=="min"
                    push!(J[temp[argmin(cliquesize[temp])]], i-1)
                else
                    push!(J[temp[argmax(cliquesize[temp])]], i-1)
                end
            else
                push!(ncc, i-1)
            end
        end
    end
    return J,ncc
end

function clique_decomp(n::Int,supp::Vector{Vector{UInt16}};alg="MD",minimize=false)
    if alg==false
        cliques=[UInt16[i for i=1:n]]
        cql=1
        cliquesize=[n]
    else
        G=SimpleGraph(n)
        for j = 1:length(supp)
            add_clique!(G, unique(supp[j]))
        end
        if alg=="NC"
            cliques,cql,cliquesize=max_cliques(G)
        else
            cliques,cql,cliquesize=chordal_cliques!(G, method=alg, minimize=minimize)
        end
    end
    uc=unique(cliquesize)
    sizes=[sum(cliquesize.== i) for i in uc]
    println("------------------------------------------------------")
    println("The clique sizes of varibles:\n$uc\n$sizes")
    println("------------------------------------------------------")
    return cliques,cql,cliquesize
end

function clique_decomp(n::Int,m::Int,d::Int,dg::Vector{Int},supp::Vector{Vector{Vector{UInt16}}};alg="MD",minimize=false)
    if alg==false
        cliques=[UInt16[i for i=1:n]]
        cql=1
        cliquesize=[n]
    else
        G=SimpleGraph(n)
        for i=1:m+1
            if i==1||d==ceil(Int, dg[i-1]/2)
                for j = 1:length(supp[i])
                    add_clique!(G, unique(supp[i][j]))
                end
            else
                temp=copy(supp[i][1])
                for j=2:length(supp[i])
                    append!(temp, supp[i][j])
                end
                add_clique!(G, unique(temp))
            end
        end
        if alg=="NC"
            cliques,cql,cliquesize=max_cliques(G)
        else
            cliques,cql,cliquesize=chordal_cliques!(G, method=alg, minimize=minimize)
        end
    end
    uc=unique(cliquesize)
    sizes=[sum(cliquesize.== i) for i in uc]
    println("------------------------------------------------------")
    println("The clique sizes of varibles:\n$uc\n$sizes")
    println("------------------------------------------------------")
    return cliques,cql,cliquesize
end
