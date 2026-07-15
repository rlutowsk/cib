#
# CIB: Cofinite Integral Braces in GAP
#
# Implementations
#

# first the cryst extensions
InstallMethod(IsIntegralAffineCrystGroup,
	"for affine crystallographic groups",
	[IsAffineCrystGroupOnLeftOrRight],
	function(agrp)
		return IsStandardAffineCrystGroup(agrp) and IsIntegerMatrixGroup( PointGroup(agrp) );
	end
);

BindGlobal("CIBVectorSystemContextFamily",
           NewFamily("CIBVectorSystemContextFamily"));


BindGlobal("CIBVectorSystemContextType",
           NewType(CIBVectorSystemContextFamily,
                   IsCIBVectorSystemContext) );

InstallMethod(ViewString,
	"for CIB vector system contexts",
	[IsCIBVectorSystemContext],
	function(obj)
		return Concatenation(
			"<CIB Vector System Context on ",
			String(UnderlyingGroup(obj)),
			">"
		);
	end
);

BindGlobal("CIBVectorSystemContextDataFamily",
           NewFamily("CIBVectorSystemContextDataFamily"));


BindGlobal("CIBVectorSystemContextDataType",
           NewType(CIBVectorSystemContextDataFamily,
                   IsCIBVectorSystemContextData and IsDataObjectRep) );

InstallMethod(ViewString,
	"for CIB vector system contexts raw data",
	[IsCIBVectorSystemContextData],
	function(obj)
		return "<CIB Vector System Context Data>";
	end
);

InstallGlobalFunction(CIBVectorSystemContext,
	function(grp, exp)
		local ctx, raw;

		ctx := Objectify( CIBVectorSystemContextType, rec() );

		raw := CIBVectorSystemContextCreate( exp, VectorSystem(grp)*exp, CoboundaryBasisInt(grp) );

		# if raw = fail then
        #     return fail;
        # fi;

		ctx!.data := raw;

		SetUnderlyingGroup(ctx, grp);

		SetExponent(ctx, exp);

		return ctx;
	end
);

InstallGlobalFunction(ClearCIBMaxExponentCache,
	function()
		CIBMaxExponentCache.abelian_groups := [];
		CIBMaxExponentCache.exponents := [];
	end
);

CoboundaryBasisIntGensOnRight := function( generators )
    local equation, deg, id, long, snf, len, res;

    deg := Size( generators[1] );
    id  := IdentityMat( deg );
    # construct the equation
    # the result of this should be a matrix of the form
    # [ x_1 - 1, x_2 - 1, ..., x_n - 1 ]
    equation := TransposedMat(
				Concatenation(
					List(
						generators,
						x->TransposedMat(x-id)
					)
				)
			);
    snf := SmithNormalFormIntegerMatTransforms( equation );
    long := (snf.coltrans^-1){[1..snf.rank]};
    len := Length( generators )-1;
    res := List(long, row->List([0..len], i->row{[i*deg+1..(i+1)*deg]}));

    return res;
end;

CoboundaryBasisIntOnRightByGens := function( generators )
    local gen_vec, gen_grp, gen_len, orbit_grp, orbit_vec, i, j, n, dd;

    gen_vec := TransposedMat( CoboundaryBasisIntGensOnRight( generators ) );
    gen_grp := List( generators, MutableCopyMat );
    gen_len := Length( gen_grp );

    orbit_grp := [ One(gen_grp[1]) ];
    orbit_vec := [ 0*gen_vec[1] ];

    i := 1;
    while i<=Length( orbit_grp ) do
        dd := orbit_grp[i];
        for j in [1..gen_len] do
            n := dd * gen_grp[j];
            if not n in orbit_grp then
                Add( orbit_grp, n );
                Add( orbit_vec, orbit_vec[i] * gen_grp[j] + gen_vec[j]);
            fi;
        od;
        i := i+1;
    od;
    # sorting is important:
	# the elements of the result are in the same order as
	# sorted elements of the point group
    SortParallel( orbit_grp, orbit_vec );

    return TransposedMat( orbit_vec );
end;


InstallMethod(CoboundaryBasisInt,
    "for affine crystallographic groups acting on right",
	[IsIntegralAffineCrystGroup],
    function(agrp)
        local res;
        if IsAffineCrystGroupOnRight( agrp ) then
            res := CoboundaryBasisIntOnRightByGens( GeneratorsOfGroup( PointGroup( agrp ) ) );
        else
            res := CoboundaryBasisIntOnRightByGens( List( GeneratorsOfGroup( PointGroup( agrp ) ), TransposedMat ) );
        fi;
		return Immutable(res);
    end
);

InstallGlobalFunction("AllAbelianPGroups",
	function( size )
		local p, n, ps;

		if not IsPrimePowerInt( size ) then
			return fail;
		fi;
		p  := PrimeDivisors( size )[1];
		n  := LogInt( size, p );
		ps := Partitions( n );
		return List(ps,
			part -> AbelianGroup( List(part, i->p^i) )
		);
	end
);

InstallGlobalFunction("AllAbelianGroups",
	function( size )
		local ppi, abp;

		if not IsPosInt(size) then
			return fail;
		fi;
		if size = 1 then
			return [ TrivialGroup() ];
		fi;
		ppi := PrimePowersInt( size );
		ppi := List([1..Size(ppi)/2], i->ppi[2*i-1]^ppi[2*i]);
		abp := List(ppi, pp->AllAbelianPGroups( pp ));
		return List(Cartesian(abp), DirectProduct);
	end
);

CIB.GetAbelianSubgroupTypes := function(invariants)
    local p_parts, p, inv, p_types, all_p_types, current_combos, new_combos, type, c;

    # 1. Rozbijamy niezmienniki na składowe p-prymarne
    # Przykład: [2, 2, 9, 64] -> parts[2] = [2, 2, 64], parts[3] = [9]
    p_parts := rec();
    for inv in invariants do
        for p in PrimeDivisors(inv) do
            if not IsBound(p_parts.(String(p))) then
                p_parts.(String(p)) := [];
            fi;
            Add(p_parts.(String(p)), p^Valuation(inv, p));
        od;
    od;

    all_p_types := [];

    # 2. Dla każdej części p-prymarnej generujemy możliwe pod-niezmienniki
    for p in RecNames(p_parts) do
        inv := p_parts.(p);
        # Generujemy wszystkie ciągi (mu_1, ..., mu_k) takie, że 0 <= mu_i <= e_i
        # Gdzie e_i to wykładniki p^e_i
        p_types := Set(
			List(
				Cartesian(
					List(inv, x -> [0..LogInt(x, Int(p))])
				),
            	vals -> SortedList(List(Filtered(vals, v -> v > 0), v -> Int(p)^v))
			)
		);
        Add(all_p_types, p_types);
    od;

    # 3. Składamy części p-prymarne (Produkt kartezjański typów)
    current_combos := [[]];
    for p_types in all_p_types do
        new_combos := [];
        for c in current_combos do
            for type in p_types do
                Add(new_combos, Concatenation(c, type));
            od;
        od;
        current_combos := new_combos;
    od;

    # 4. Sortujemy każdy wynik i zwracamy unikalne zestawy
    return Set(List(current_combos, x -> SortedList(x)));
end;

IsIsomorphicSubgroup := function (G, H)
    local p, r, Q, h, target_size, iso_it;

    target_size := Size(H);

    # necessary condition (orders)
    if Size(G) mod target_size <> 0 then return false; fi;

    # p-subgroup case
    p := IsPrimePowerInt(target_size);
    if p <> false then
        Info(InfoCIB, 2, "H is a p-group. Searching in Sylow subgroup.");
        return IsomorphicSubgroups(SylowSubgroup(G, p), H : findall:=false) <> [];
    fi;

    # 3. Wykorzystanie radykału (Twój pomysł, ale szybciej)
    r := SolvableRadical(G);
    if target_size <= Size(r) then
        if IsomorphicSubgroups(r, H : findall:=false) <> [] then
            return true;
        fi;
    fi;

    # 4. Zamiast ręcznego liftingu, użyj Iteratora (Oszczędność pamięci)
    # IsomorphicSubgroups dla grup nierozwiązalnych używa techniki
    # "lifting through the solvable radical" automatycznie.
    Info(InfoCIB, 2, "Starting general search via IsomorphicSubgroups iterator.");

    # IsomorphicSubgroups(G, H) zwraca listę, ale możemy użyć wersji
    # która nie wylicza wszystkiego na raz (jeśli dostępna w Twojej wersji GAP)
    if IsomorphicSubgroups(G, H : findall:=false) <> [] then
        return true;
    fi;

    return false;
end;

MaximalSolvableSubgroupClassReps := function( grp )
	local rad, hom, img;

	if IsSolvableGroup( grp ) then
		return [ grp ];
	fi;

	rad := SolvableRadical( grp );
	hom := NaturalHomomorphismByNormalSubgroup( grp, rad );
	return List( MaximalSolvableSubgroups(Image( hom )), x->PreImage( hom, x ) );
end;

CIB.BraceMaxExponents := function(grp)
	local subgroups, quotients, allAbelianGroups, exp, agrp, availableIds, aut, sg, qg, allowedIds, exps, i, k, img, max, sub, flt;

	if not IsSolvableGroup(grp) then
		return [];
	fi;

	subgroups := [];
	quotients := [];
	for sg in Filtered( NormalSubgroups(grp), IsAbelian ) do
		qg := grp/sg;
		if ForAny( quotients, q->IsomorphismGroups(q, qg)<>fail ) then
			continue;
		fi;
		Add( subgroups, sg );
		# Add( quotients, Group( GeneratorsSmallest( qg ) ) );
		Add( quotients, Image( SmallerDegreePermutationRepresentation( Image( IsomorphismPermGroup( qg )  ) ) ) );
	od;

	allAbelianGroups := AllAbelianGroups( Size(grp) );

	exps := [];
	for exp in Reversed( DivisorsInt( Size(grp) ) ) do
		Info( InfoCIB, 3, "Exponents' list: ", exps);
		if ForAny(exps, e->IsInt(e/exp)) then
			continue;
		fi;
		Info(InfoCIB, 3, "Checking for exponent: ", exp);
		for agrp in Filtered( allAbelianGroups, g -> Exponent(g) = exp ) do
			if exp in exps then
				break;
			fi;
			Info(InfoCIB, 3, "Checking abelian group A: ", AbelianInvariants(agrp));

			availableIds := CIB.GetAbelianSubgroupTypes( AbelianInvariants(agrp) );
			aut := AutomorphismGroup( agrp );
			img := Image( SmallerDegreePermutationRepresentation( NiceObject( aut ) ) );
			Info(InfoCIB, 3, "Aut(A) is of order: ", Size(aut));

			max := List( MaximalSolvableSubgroupClassReps( img ), x->Image( IsomorphismPcGroup( x ) ) );
			Sort(max, function(x,y) return Size(x)< Size(y); end );
			Info(InfoCIB, 3, "Maximal solvable subgroups of Aut(A) have orders: ", List(max, Size));

			# sub := List( max, x -> SubgroupsSolvableGroup( x, rec( consider := SizeConsiderFunction( Size(grp) ) ) ) );
			# Info(InfoCIB, 3, "Subgroups to consider in each maximal subgroup: ", List(sub, Size));
			sub := [];

			for i in [1..Size(max)] do
				sub[i] := SubgroupsSolvableGroup( max[i], rec( consider := SizeConsiderFunction( Size(grp) ) ) );
				Info( InfoCIB, 3, "Calculated subgroups to consider of maximal subgroup number ", i);

				for k in [1..Size(subgroups)] do
					sg := subgroups[k];
					Info(InfoCIB, 3, "Checking subgroup number ", k,", : ", AbelianInvariants(sg));
					if not AbelianInvariants(sg) in availableIds then
						continue;
					fi;

					# Place to put some filters which are cheap and may be applied before the isomorphism test.
					# The simple examples below (number 1, 2 and 3 - in more detail) are already done in the `IsomorphismGroups',
					# but for some groups 2 and 3 are the last step of the calculations, so we try them here first.
					#
					# We know that the groups involved are solvable.

					# 1. Filter by the order of the group
					flt := Filtered( sub[i], s->Size(s) = Size(quotients[k]) );
					Info(InfoCIB, 4, "Passed size filter: ", Size(flt));

					# 2. Filter by abelianization - done by `IsomorphismGroups'
					flt := Filtered( flt, s->AbelianInvariants(s) = AbelianInvariants(quotients[k]) );
					Info(InfoCIB, 4, "Passed abelianization filter: ", Size(flt));

					# 3. Filter by number of conjugacy classes - done by `IsomorphismGroups'
					flt := Filtered( flt, s->NrConjugacyClasses(s) = NrConjugacyClasses(quotients[k]) );
					Info(InfoCIB, 4, "Passed number of conjugacy classes filter: ", Size(flt));

					# n. Filter by isomorphism (the most expensive one, so we do it at the end)
					if First( flt, s->IsomorphismGroups( quotients[k], s ) <> fail ) <> fail then
						Add(exps, exp);
						break;
					fi;
				od;
			od;
		od;
	od;
	return exps;
end;

CIB.CacheInit := function( grp )
	local size;

	size := Size(grp);

	if not IsBound(CIB.cache) then
		CIB.cache := rec( abelian_groups := [], exponents := [] );
	else
		if not IsBound(CIB.cache.exponents) then
			CIB.cache.exponents := [];
		fi;
		if not IsBound(CIB.cache.abelian_groups) then
			CIB.cache.abelian_groups := [];
		fi;
	fi;
	if not IsBound(CIB.cache.exponents[ size ]) then
		CIB.cache.exponents[ size ] := [];
	fi;
	if not IsBound(CIB.cache.abelian_groups[ size ]) then
		CIB.cache.abelian_groups[ size ] := AllSmallGroups( size, IsAbelian );
	fi;
	return CIB.cache;
end;

CIB.CacheTry := function( grp )
	local id;

	id := IdGroup(grp);

	if IsBound( CIB.cache.exponents[id[1]][id[2]] ) then
		return CIB.cache.exponents[id[1]][id[2]];
	fi;
	return fail;
end;

CIB.CacheGetAbelianGroups := function( grp )
	return CIB.cache.abelian_groups[Size(grp)];
end;

CIB.CacheSetExponents := function( grp, exps )
	local id;

	id := IdGroup(grp);

	CIB.cache.exponents[id[1],id[2]] := exps;
end;

CIB.BraceMaxExponentsSmallGroup := function(grp)
	local exp, allAbelianGroups, agrp, subgroups, availableIds, aut, sg, allowedIds, cache, exps;

	if not IsSolvableGroup(grp) then
		return [];
	fi;

	subgroups := Filtered( NormalSubgroups(grp), IsAbelian );

	CIB.CacheInit( grp );

	cache := CIB.CacheTry( grp );

	if cache<> fail then
		return cache;
	fi;

	allAbelianGroups := CIB.CacheGetAbelianGroups( grp );

	exps := [];

	for exp in Reversed( DivisorsInt( Size(grp) ) ) do
		Info( InfoCIB, 3, "Exponents' list: ", exps);
		if ForAny(exps, e->IsInt(e/exp)) then
			continue;
		fi;
		Info(InfoCIB, 3, "Checking for exponent: ", exp);
		for agrp in Filtered( allAbelianGroups, g -> Exponent(g) = exp ) do
			if exp in exps then
				break;
			fi;
			Info(InfoCIB, 3, "Checking abelian group with id: ", IdGroup(agrp));
			availableIds := Unique(List(
				# SubgroupsSolvableGroup(agrp),
				List(ConjugacyClassesSubgroups( agrp ), Representative),
				IdGroup
			));
			# Info(InfoCIB, 3, "Possible ids: ", availableIds);
			for sg in subgroups do
				Info(InfoCIB, 3, "Checking subgroup: ", IdGroup(sg));
				if not IdGroup(sg) in availableIds then
					continue;
				fi;
				# Display(IdGroup(grp/sg));
				aut := AutomorphismGroup( agrp );
				# check the obvious condition on the order of the group
				if not IsInt(Size(aut)*Size(sg)/Size(grp)) then
					continue;
				fi;
				# if it is too large, then return exp as possible answer
				if Size( aut ) > 10^4 then
					Add(exps, exp);
					break;
				fi;
				allowedIds := SSortedList(
					Filtered(
						List(
							ConjugacyClassesSubgroups(
								NiceObject( aut )
							),
							Representative
						),
						cs->Size(cs)*Size(sg) = Size(grp)
					),
					IdGroup
				);
				Info(InfoCIB, 3, "Subgroups of Aut(grp): ", allowedIds);
				if IdGroup(grp/sg) in allowedIds then
					Add(exps, exp);
					break;
				fi;
			od;
		od;
	od;
	CIB.CacheSetExponents( grp, exps );
	return Immutable(exps);
end;

InstallMethod(BraceMaxExponents,
	"find maximal possible exponents of additive groups of braces with given multiplicative one",
	[IsGroup and IsFinite],
	function(grp)
		local img;
		if not IsSolvableGroup(grp) then
			return [];
		fi;
		if IsPcGroup( grp ) then
			img := grp;
		else
			img := Image( IsomorphismPcGroup( grp ) );
		fi;
		if IdGroupsAvailable( Size( img ) ) then
			return CIB.BraceMaxExponentsSmallGroup( img );
		fi;
		return CIB.BraceMaxExponents( img );
	end
);

InstallMethod(BraceMaxExponents,
	"find maximal possible exponents of additive groups of braces with multiplicative one given by point group of a crystallographic group",
	[IsAffineCrystGroupOnLeftOrRight],
	function(sgrp)
		return BraceMaxExponents( PointGroup(sgrp) );
	end
);

CIB.GetZmodnZGroup := function(list, m)
    local l, s, i, j;

	l := list mod m;
	if not 0*l[1] in l then
        return fail;
    fi;
    s := Size(l);
    for i in [1..s] do
        if not (-l[i] mod m) in l then
            return fail;
        fi;
        if not (2*l[i] mod m) in l then
            return fail;
        fi;
        for j in [i+1..s] do
            # check uniqueness of elements of the list
            if l[i] = l[j] then
                return fail;
            fi;
            # check addition
            if not ((l[i]+l[j]) mod m) in l then
                return fail;
            fi;
        od;
    od;
    return l;
end;

InstallMethod(IsZmodnZGroup,
	"for list of vectors and modulus",
	[IsMatrix, IsPosInt],
	function(list, m)
		return CIB.GetZmodnZGroup(list, m) <> fail;
	end
);

CIB.CofiniteIntegralBraceVectorSystemsGAP := function( agrp )
	local c, d, g, l, b, i, j, k, cp, vs, cnt, pos, rank, pgrp, cib, exp, curr_l, l_mod, res, den;

	b := CoboundaryBasisInt( agrp );
	d := DegreeOfMatrixGroup( agrp ) - 1;

	exp := BraceMaxExponents( agrp );

	pgrp := PointGroup( agrp );

	pos := Position( Set(pgrp), One(pgrp) );
	if ForAny(b, v->ForAny(v[pos], x->x<>0)) then
		Error("wrong assumptions - identity not mapped to zero");
	fi;

	rank  := Size(b);
	cib := [];

	den := Lcm( List(Flat(VectorSystem(agrp)), DenominatorRat) );

	for g in exp do
		g  := Lcm( g, den );
		vs := VectorSystem(agrp) * g;

		# Preliminary check
		# if ForAny( Flat(vs), x->not IsInt(x) ) then
		# 	continue;
		# fi;

		Info(InfoCIB, 2, g^rank, " items to process");

		# curr_l holds the current state (without modulo)
		curr_l := ShallowCopy(vs);

		# coordinate vector
		c := ListWithIdenticalEntries(rank, 0);

		while true do

			res := CIB.GetZmodnZGroup( curr_l, g );
			if res<>fail then
				Info(InfoCIB, 3, "Found: ", c, curr_l, res);
				Add(cib, res/g);
			fi;

			# here we change the state (coordinate vector)
			pos := rank;

			# We search for a position to increment (from right to left)
			# If a digit reached max (g-1), we reset it and carry over further
			while pos > 0 and c[pos] = g - 1 do
				# RESET POSITION:
				# Since c[pos] was g-1, we added (g-1) times vector b[pos] to the sum.
				# We must now subtract it to return to the state as if c[pos] was 0.
				curr_l := curr_l - (g-1) * b[pos];

				c[pos] := 0;
				pos := pos - 1;
			od;

			# If we've processed all positions, we're done with this exponent
			if pos = 0 then
				break;
			fi;

			# INCREMENT:
			# We increment the counter at the found position and update the vector sum
			c[pos] := c[pos] + 1;
			curr_l := curr_l + b[pos];
		od;
	od;
	return Set( cib );
end;

CIB.FixDenominators := function( agrp, g )
	local gens, ind, P, rhs, lhs, sol, sgens, i, tmp, d, d1, nd, vec, sgrp;

	P   := PointGroup( agrp );
	gens:= GeneratorsOfGroup(P);
	ind := List( gens, x->Position(SSortedList(P),x) );
	rhs := Concatenation( VectorSystem(agrp){ind} );
	lhs := List( CoboundaryBasisInt(agrp), row->Concatenation(row{ind}) );
	sol := SolveInhomEquationsModZ(lhs, g*rhs, true);
	if sol[1] = [] then
		return fail;
	fi;
	lhs := List( CoboundaryBasisInt(agrp), row->row{ind} );
	rhs := VectorSystem(agrp){ind} - Sum([1..Size(sol[1,1])], i->sol[1,1][i]*lhs[i])/g;

	sgens := [];
	d := DegreeOfMatrixGroup(P);
	d1:= d+1;
	nd:= [1..d];
	for i in [1..Size(gens)] do
		tmp := NullMat(d1,d1);
		tmp{nd}{nd} := gens[i];
		tmp[d1]{nd} := List(rhs[i], FractionModOne);
		tmp[d1,d1]  := 1;
		Add( sgens, tmp );
	od;
	for vec in TranslationBasis( agrp ) do
		tmp := IdentityMat(d1);
		tmp[d1]{nd} := vec;
		Add( sgens, tmp );
	od;
	sgrp := AffineCrystGroupOnRight( sgens );
	SetBraceMaxExponents( sgrp, BraceMaxExponents( agrp ) );
	return sgrp;
end;

CIB.CofiniteIntegralBraceVectorSystemsByContext := function(agrp)
	local exp, cib, res, g, name, ctx, rank, den, good;

	exp := BraceMaxExponents( agrp );

	cib := [];

	if HasName(agrp) then
		name := Concatenation(Name(agrp), ": ");
	else
		name := "";
	fi;

	rank := Size(CoboundaryBasisInt(agrp));
	den := Lcm( List(Flat(VectorSystem(agrp)), DenominatorRat) );

	for g in exp do
		if g <> Lcm(g, den) then
			good := CIB.FixDenominators( agrp, g );
		else
			good := agrp;
		fi;
		if good = fail then
			continue;
		fi;
		Info(InfoCIB, 3, name, g^rank, " items to process");
		ctx := CIBVectorSystemContext( good, g );
		if ctx = fail then
			continue;
		fi;
		res := CofiniteIntegralBraceVectorSystems( ctx );
		Info(InfoCIB, 3, name, "exponent ", g, ": found ", Size(res), " vector systems.");
		Append(cib, res );
	od;
	# SetCofiniteIntegralBraceVectorSystems( agrp, cib );
	return Set( cib );
end;

CIB.BySubgroupPointGroup := function( sgrp, spg, cnt... )
    local cb, ind, ssg, cib, vs, c, sol, cand, res, rk, m, solve, one, basis, coeffs, x;

    ind := List( SSortedList(spg), g -> Position( SSortedList( PointGroup(sgrp) ), g ) );
    ssg := PreImage( PointHomomorphism(sgrp), spg );
	if PointGroup(ssg) <> spg then
        Error("The subgroup is not a subgroup of the point group of the group. This should not happen.");
    fi;

    if InfoLevel(InfoCIB) > 2 and Size(cnt)>0 and IsInt(cnt[1]) then
        SetName(ssg, StringFormatted("{} / {}", Name(sgrp), cnt[1] ));
        Info( InfoCIB, 3, Name(ssg), ": starting calculations");
    fi;

    vs  := VectorSystem( sgrp );
    cib := List( CIB.CofiniteIntegralBraceVectorSystemsByContext( ssg ), x->Concatenation(x-vs{ind}));

    cb  := List( CoboundaryBasisInt( sgrp ), x -> Concatenation( x{ind} ) );
    sol := [];
    m   := Size( PointGroup(sgrp) );
    for c in cib do
		solve := SolveInhomEquationsModZ(cb, c, true);
		one   := solve[1] * m;
		basis := solve[2];
		if basis = [] then
			Append( sol, one );
			continue;
		fi;
		# taken from NullspaceModN code
		coeffs:= Cartesian( ListWithIdenticalEntries( Length( basis ), [ 0 .. m - 1 ] ) );
		for x in one do
        	Append( sol, Set(coeffs, cf->(x + cf * basis) mod m ) );
		od;
    od;
    res := [];
    cb  := CoboundaryBasisInt( sgrp );
    rk  := Size( cb );
    vs  := vs * m;
    for c in sol do
        cand := CIB.GetZmodnZGroup( vs + Sum([1..rk], i->c[i]*cb[i]), m );
        if cand<>fail then
            Add( res, cand/m );
        fi;
    od;
	# Error("CHECK");
    return res;
end;

InstallMethod(CofiniteIntegralBraceVectorSystems,
    "for affine crystallographic groups acting on right",
    [IsAffineCrystGroupOnRight and IsStandardAffineCrystGroup],
    function(agrp)
		local syl, res, p, pd, d, sg, iso, img;

		if not IsSolvableGroup( agrp ) then
			return [];
		fi;

		pd := PrimeDivisors( Size( PointGroup( agrp ) ) );
		if pd = [] then
			return [[ ListWithIdenticalEntries( DegreeOfMatrixGroup(agrp)-1, 0 ) ]];
		fi;

		res := [];
		d   := 0;
		if Size(pd) = 1 then
			# handle p-groups
			p := pd[1];
			Info( InfoCIB, 1, Name(agrp), ": calculating CIBs for ", p, "-group of order ", Size(PointGroup(agrp)) );
			if IsOddInt( p ) or Size(PointGroup(agrp)) <= 32 then
				return CIB.CofiniteIntegralBraceVectorSystemsByContext( agrp );
			fi;
			iso := IsomorphismPcGroup( PointGroup(agrp) );
			img := Image( iso );
			for sg in Concatenation(
				List(
					Filtered(
						ConjugacyClassesSubgroups( img ),
						x->Size(Representative(x))=16
					),
					List)
				) do
				d := d+1;
				Append( res, CIB.BySubgroupPointGroup( agrp, PreImage(iso, sg), d ) );
				res := Unique( res );
			od;
			return res;
		fi;
		# take first odd prime and calculate with this one
		p := pd[2];
		Info( InfoCIB, 1, Name(agrp), ": calculating CIBs using ", p, "-Sylow subgroups of order ", Size(SylowSubgroup(PointGroup(agrp), p)) );

		for syl in SylowSubgroup( PointGroup( agrp ), p )^PointGroup( agrp ) do
			d := d+1;
			Append( res, CIB.BySubgroupPointGroup( agrp, syl, d ) );
			res := Unique( res );
		od;

		return res;
    end
);

CIB.ExtendedVectorSystems := function( sgrp, ssg, spg, cibs )
	local cb, ind, vs, c, sol, cib, cand, res, rk, m, solve, one, basis, coeffs, x, cb_full, ctx;

	ind := List( SSortedList(spg), g -> Position( SSortedList( PointGroup(sgrp) ), g ) );

	# if PointGroup(ssg) <> spg then
	# 	Error("The subgroup is not a subgroup of the point group of the group. This should not happen.");
	# fi;

	# if InfoLevel(InfoCIB) > 2 and Size(cnt)>0 and IsInt(cnt[1]) then
	# 	SetName(ssg, StringFormatted("{} / {}", Name(sgrp), cnt[1] ));
	# 	Info( InfoCIB, 3, Name(ssg), ": starting calculations");
	# fi;
	Info( InfoCIB, 3, "extending vector systems");

	vs  := VectorSystem( sgrp );
	cib := List( cibs, x->Concatenation(x-vs{ind}));

	cb  := List( CoboundaryBasisInt( sgrp ), x -> Concatenation( x{ind} ) );
	sol := [];
	
    cb_full := CoboundaryBasisInt( sgrp );
    rk  := Size( cb );
    m   := Size( PointGroup(sgrp) );
    vs  := vs * m;
    res := [];
    Info( InfoCIB, 3, "will run ", Size(cib), " times SolveInhomEquationsModZ" );
	for c in cib do
		solve := SolveInhomEquationsModZ(cb, c, true);
        Info( InfoCIB, 4, "solved");
		one   := Filtered(solve[1] * m, x->ForAll(x,IsInt) );
		basis := solve[2];
		if basis = [] then
			# CHECK IT
			# Idea as follows: if one is not an integer vector, skip it
			# Append( sol, one );
            for c in one do
    		    cand := CIB.GetZmodnZGroup( vs + Sum([1..rk], i->c[i]*cb_full[i]), m );
	    	    if cand<>fail then
		    	    Add( res, cand/m );
    		    fi;
	        od;
			continue;
		fi;
		# taken from NullspaceModN code
		#coeffs:= Cartesian( ListWithIdenticalEntries( Length( basis ), [ 0 .. m - 1 ] ) );
		for x in one do
            #for c in Set(coeffs, cf->(x + cf * basis) mod m ) do
            #   cand := CIB.GetZmodnZGroup( vs + Sum([1..rk], i->c[i]*cb_full[i]), m );
	    	#   if cand<>fail then
		    # 	    Add( res, cand/m );
    	    #	fi;
            #od;
            # try with context
            ctx := Objectify( CIBVectorSystemContextType, rec() );
            #Error("check");
            ctx!.data := CIBVectorSystemContextCreate( m, (vs + Sum([1..rk], i->x[i]*cb_full[i])) mod m, List(basis, b->Sum([1..rk],i->b[i]*cb_full[i])) );
            SetUnderlyingGroup( ctx, sgrp );
            SetExponent( ctx, m );
            Append( res, CofiniteIntegralBraceVectorSystems( ctx ) );
            Info( InfoCIB, 5, "checked" );
		od;
	od;
	#cb  := CoboundaryBasisInt( sgrp );
	#rk  := Size( cb );
	#vs  := vs * m;
    #Display(res);	
	return res;
end;

CIB.ByPSubgroupPointPGroup := function( sgrp, depth )
    local pg, ord, spg, ssg, cib, res;

	pg  := PointGroup( sgrp );
	if not IsPGroup( pg ) or Size(pg) = 1 then
		Error("Point group of <sgrp> must be a (non-trivial) p-group");
	fi;
	ord := Size( pg );
	if ord <= 32 or Length(FactorsInt(ord))<=2 then
		Info(InfoCIB, 3, "calculating CIBs directly");
		Info(InfoCIB, 4, VectorSystem(sgrp) );
		return CIB.CofiniteIntegralBraceVectorSystemsByContext( sgrp );
	fi;

	res := [];
	Info( InfoCIB, 2, depth, ": Checking ", Size(MaximalNormalSubgroups(pg)), " normal subgroups");
	for spg in MaximalNormalSubgroups( pg ) do
		ssg := PreImage( PointHomomorphism(sgrp), spg );
		cib := CIB.ByPSubgroupPointPGroup( ssg, depth+1 );
		Append( res, CIB.ExtendedVectorSystems( sgrp, ssg, spg, cib ) );
	od;
	# Error("CHECK");
    return SSortedList( res );
end;

CIB.ByPSubgroupPointPQGroup := function( sgrp, cnt... )
    local pg, pd, max, syl, cib, ssg, res;

	pg  := PointGroup( sgrp );
	pd  := PrimeDivisors( Size(pg) );
	max := Maximum( List(pd, p->Size(SylowSubgroup( pg, p )) ) );

	res := [];
	for syl in SylowSubgroup( pg, PrimeDivisors(max)[1] )^pg do
		ssg := PreImage( PointHomomorphism(sgrp), syl );
		cib := CIB.ByPSubgroupPointPGroup( ssg, 0 );
		# Error("CHECK");
    	Append( res, CIB.ExtendedVectorSystems( sgrp, ssg, syl, cib ) );
	od;
	return SSortedList( res );
end;

CIB.CofiniteIntegralBraceVectorSystemsPart := function(agrp)
	local pd;

	Info( InfoCIB, 1, Name(agrp), ": calculating CIBs for point group of order ", Size(PointGroup(agrp)) );

	if not IsSolvableGroup( agrp ) then
		return [];
	fi;

	pd := PrimeDivisors( Size( PointGroup( agrp ) ) );
	if pd = [] then
		return [[ ListWithIdenticalEntries( DegreeOfMatrixGroup(agrp)-1, 0 ) ]];
	fi;

	if Size(pd) = 1 then
		# handle p-groups
		return CIB.ByPSubgroupPointPGroup( agrp, 0 );
	fi;

	return CIB.ByPSubgroupPointPQGroup( agrp );
end;

InstallMethod(CofiniteIntegralBraceVectorSystems,
    "for affine crystallographic groups acting on left",
    [IsAffineCrystGroupOnLeft and IsStandardAffineCrystGroup],
	function(agrp)
		return CofiniteIntegralBraceVectorSystems( TransposedMatrixGroup(agrp) );
	end
);

InstallMethod(CofiniteIntegralBraceVectorSystems,
	"for CIB vector system contexts",
	[IsCIBVectorSystemContext],
	function(ctx)
		local res, exp;
		res := [];
		exp := Exponent(ctx);
		repeat
			# Info( InfoCIB, 4, "Checking coordinates ", CIBVectorSystemContextRawCoordinates(ctx!.data) );
			if IsGroupByCIBVectorSystemContextRaw(ctx!.data) then
				Add( res, CIBVectorSystemContextRawMat(ctx!.data, 0) / exp );
			fi;
		until not CIBVectorSystemContextRawInc(ctx!.data);
		return res;
	end
);

InstallMethod(VectorSystem,
    "for affine crystallographic groups acting on right",
    [IsAffineCrystGroupOnRight],
    function(agrp)
		local pg, ph, d, vs;

		pg := PointGroup(agrp);
		ph := PointHomomorphism(agrp);
		d  := DegreeOfMatrixGroup(pg);

		# return List( SSortedList(pg), x->List(PreImagesRepresentative( ph, x )[d+1]{[1..d]}, FractionModOne) );
        vs := (List( SSortedList(pg), x->PreImagesRepresentative( ph, x )[d+1]{[1..d]} ) * Size(pg) ) mod Size(pg);
        return vs / Size(pg);
    end
);

InstallMethod(VectorSystem,
    "for affine crystallographic groups acting on left",
    [IsAffineCrystGroupOnLeft],
	function(agrp)
		return VectorSystem( TransposedMatrixGroup(agrp) );
	end
);

InstallMethod(GeneratorsOfIntegralAffineCrystGroupOnRight,
	"for integral affine crystallographic groups acting on right",
	[IsIntegralAffineCrystGroup and IsAffineCrystGroupOnRight],
	function(agrp)
		return GeneratorsOfGroup(agrp);
	end
);

InstallMethod(GeneratorsOfIntegralAffineCrystGroupOnRight,
	"for list of matrices generating affine crystallographic groups acting on right",
	[IsList],
	function(gens)
		local d, f;

		if not ForAll(gens, IsAffineMatrixOnRight) then
			Error("Generators must be affine matrices acting on right");
		fi;

		d := Size(gens[1]);
		f := function(x)
			local m;
			m := IdentityMat(d);
			m[d,x] := 1;
			return m;
		end;
		return Concatenation( gens, List([1..d-1], f) );
	end
);

InstallMethod(GeneratorsOfIntegralAffineCrystGroupOnLeft,
	"for integral affine crystallographic groups acting on left",
	[IsIntegralAffineCrystGroup and IsAffineCrystGroupOnLeft],
	function(agrp)
		return GeneratorsOfGroup(agrp);
	end
);

InstallMethod(GeneratorsOfIntegralAffineCrystGroupOnLeft,
	"for list of matrices generating affine crystallographic groups acting on left",
	[IsList],
	function(gens)
		local d, f;

		if not ForAll(gens, IsAffineMatrixOnLeft) then
			Error("Generators must be affine matrices acting on left");
		fi;

		d := Size(gens[1]);
		f := function(x)
			local m;
			m := IdentityMat(d);
			m[x,d] := 1;
			return m;
		end;
		return Concatenation( gens, List([1..d-1], f) );
	end
);

InstallMethod(CofiniteIntegralBracesRepsGenerators,
	"for integral affine crystallographic group on right",
	[IsIntegralAffineCrystGroup and IsAffineCrystGroupOnRight],
	function(agrp)
        local act, orb, cib, N, Ngens, Npair, Nperm, P, Plist, ind, fun, d, order, nr, nc, nd, d1;

        P     := PointGroup( agrp );
		d     := DegreeOfMatrixGroup( P );
		order := Size(P);
		Plist := SSortedList(P);
		N     := NormalizerInGLnZ( P );
		Ngens := [ One(N) ];
		if HasGeneratorsSmallest(N) then
			Append(Ngens, GeneratorsSmallest(N));
		else
			Append(Ngens, GeneratorsOfGroup( N ));
		fi;
		Npair := List( Ngens, n->[n^-1, n] );
		Nperm := List( Npair, n->PermList( List(Plist, e->PositionSet( Plist, n[1] * e * n[2] )) ) );

		if HasCofiniteIntegralBraceVectorSystems( agrp ) then
			cib   := CofiniteIntegralBraceVectorSystems( agrp ) * order; # cancel all denominators
		else
			cib   := CIB.CofiniteIntegralBraceVectorSystemsPart( agrp ) * order; # cancel all denominators
		fi;
		nr := [1..order];
		nc := [1..DegreeOfMatrixGroup( P )];
        act := function( list, pos )
			local img, i, j;
			img := Permuted( list, Nperm[pos] );
			img := MutableCopyMat( img * Ngens[pos] );
			for i in nr do
				for j in nc do
					img[i,j] := img[i,j] mod order;
				od;
			od;
			return img;
		end;

		orb := Orbits( N, cib, Ngens, [1..Size(Ngens)], act);
        orb := List( orb, o->First(o, x->x in cib) );

		ind := List( GeneratorsOfGroup( P ), g->Position( Plist, g ) );
		if ind = [] then
			# the point group is trivial
			ind := [1];
		fi;

		nd := [1..d];
		d1 := d+1;
		fun := function( mat, vec )
			local res, i, j;
			res := NullMat(d1,d1);
			res{nd}{nd} := mat;
			res[d1]{nd} := vec;
			res[d1,d1]  := 1;
			return res;
		end;

		orb := List(orb, o->List( ind, i->fun( Plist[i], o[i]/order ) ) );
        return orb;
	end
);

InstallMethod(CofiniteIntegralBracesRepsGenerators,
	"for integral affine crystallographic group on left",
	[IsIntegralAffineCrystGroup and IsAffineCrystGroupOnLeft],
	function(agrp)
		local tgrp, gens;
		tgrp := TransposedMatrixGroup( agrp );
		gens := CofiniteIntegralBracesRepsGenerators( tgrp );
		return List( gens, g->List(g, TransposedMat) );
	end
);

InstallMethod(SSortedGeneratorsOfGroup,
    "for group with generators",
    [IsGroup],
    function(grp)
        return SSortedList( GeneratorsOfGroup( grp ) );
    end
);
