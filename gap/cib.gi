#
# CIB: Cofinite Integral Braces in GAP
#
# Implementations
#

# first the cryst extensions
InstallMethod(IsIntegralAffineCrystGroup,
	"for affine crystallographic groups",
	[IsStandardAffineCrystGroup],
	function(agrp)
		return IsIntegerMatrixGroup( PointGroup(agrp) );
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
                   IsCIBVectorSystemContextData) );

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
        
		if raw = fail then
            return fail;
        fi;

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

InstallMethod(CoboundaryBasisInt,
    "for affine crystallographic groups acting on right",
    #[IsAffineCrystGroupOnRight and IsStandardAffineCrystGroup],
	[IsStandardAffineCrystGroup],
    function(agrp)
        local equation, deg, id, snf, long, phom, pgrp, res;
        # retrieve the point (holonomy) group
		phom := PointHomomorphism(agrp);
		pgrp := Image(phom);
        deg  := DegreeOfMatrixGroup(pgrp);
        id   := IdentityMat(deg);
		# construct the equation
        # the result of this should be a matrix of the form
        # [ x_1 - 1, x_2 - 1, ..., x_n - 1 ]
		# sorting is important:
		# the elements of the result are in the same order as
		# sorted elements of the point group
		if IsAffineCrystGroupOnRight(agrp) then
			equation := TransposedMat(
				Concatenation(
					List(
						SSortedList(pgrp),
						x->TransposedMat(x-id)
					)
				)
			);
		else
			equation := TransposedMat(
				Concatenation(
					List(
						SSortedList(pgrp),
						x->x-id
					)
				)
			);
		fi;
        # solve it
        snf := SmithNormalFormIntegerMatTransforms(equation);
        # take the first snf.rank rows of snf.coltrans^-1
        long := (snf.coltrans^-1){[1..snf.rank]};
        # split solution into a list of lists
        res := List(long, row->List([0..Size(pgrp)-1], i->row{[i*deg+1..(i+1)*deg]}));
		return Immutable(res);
    end
);

InstallMethod(BraceMaxExponents, 
	"find maximal possible exponents of additive groups of braces with given multiplicative one",
	[IsGroup and IsFinite],
	function(grp)
		local exp, allAbelianGroups, agrp, subgroups, availableIds, aut, sg, allowedIds, cache, exps;

		if not IsSolvableGroup(grp) then
			return [];
		fi;

		subgroups := Filtered( NormalSubgroups(grp), IsAbelian );

		cache:= CIBMaxExponentCache;

		if not IsBound(cache.exponents) then
			cache.exponents := [];
		fi;
		if not IsBound(cache.abelian_groups) then
			cache.abelian_groups := [];
		fi;
		if not IsBound(cache.exponents[Size(grp)]) then
			cache.exponents[Size(grp)] := [];
		fi;
		if IsBound( cache.exponents[Size(grp)][IdGroup(grp)[2]] ) then
			return cache.exponents[Size(grp)][IdGroup(grp)[2]];
		fi;
		if not IsBound(cache.abelian_groups[Size(grp)]) then
			cache.abelian_groups[Size(grp)] := AllSmallGroups( Size(grp), IsAbelian );
		fi;
		allAbelianGroups := cache.abelian_groups[Size(grp)];

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
		cache.exponents[Size(grp)][IdGroup(grp)[2]] := exps;
		return Immutable(exps);
	end
);

InstallMethod(BraceMaxExponents, 
	"find maximal possible exponents of additive groups of braces with multiplicative one given by point group of a crystallographic group",
	[IsAffineCrystGroupOnLeftOrRight],
	function(sgrp)
		return BraceMaxExponents( PointGroup(sgrp) );
	end
);

GetZmodnZGroup := function(list, m)
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

IsZmodnZGroup	:= function(list, m)
	local res;
	
	res := GetZmodnZGroup(list, m);
	return res <> fail;
end;

InstallMethod(CofiniteIntegralBraceVectorSystems,
    "for affine crystallographic groups acting on right",
    [IsAffineCrystGroupOnRight and IsStandardAffineCrystGroup],
    function(agrp)
		local c, d, g, l, b, i, j, k, cp, vs, cnt, pos, rank, pgrp, cib, exp, curr_l, l_mod, res;

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

		for g in exp do
			# Info(InfoCIB, 2, "Checking exponent ", g);
			# cp    := IteratorOfCartesianProduct(List([1..rank], i->[0..g-1]));
			# vs    := VectorSystem(agrp) * g;
			# if ForAny( Flat(vs), x->not IsInt(x) ) then
			# 	continue;
			# fi;
			# Info(InfoCIB, 2, g^rank, " items to process");
			# cnt := 0;
			# for c in cp do
			# 	cnt := cnt+1;
			# 	l:= ShallowCopy(vs);
			# 	for i in [1..rank] do
			# 		l := l + c[i]*b[i];
			# 	od;
			# 	res := GetZmodnZGroup( l, g );
			# 	if res<>fail then
			# 		Info(InfoCIB, 3, "Found: ", c);
			# 		Add(cib, res/g);
			# 	fi;
			# od;
			# continue;
			# Pre-kalkulacja wierszy dla szybszego dostępu, jeśli b jest macierzą
            # Zakładam, że b[i] zwraca i-ty wektor (wiersz), co jest standardem w GAP dla list list
            # b_rows := b; 

            vs    := VectorSystem(agrp) * g;
            
            # Wstępne sprawdzenie z Twojego kodu
            if ForAny( Flat(vs), x->not IsInt(x) ) then
                continue;
            fi;

            Info(InfoCIB, 2, g^rank, " items to process");
            
            # Inicjalizacja stanu
            # curr_l trzyma aktualną sumę liniową (bez modulo)
            curr_l := ShallowCopy(vs);
            
            # Wektor współczynników c (nasz "licznik"), startujemy od samych zer
            c := ListWithIdenticalEntries(rank, 0);
            
            while true do
                # 1. Logika sprawdzająca (Check)
                # Wykonujemy modulo tylko tutaj, na potrzeby testu, nie psując ciągłości sumy w curr_l
                # l_mod := curr_l mod g;
                
                res := GetZmodnZGroup( curr_l, g );
                if res<>fail then
					Info(InfoCIB, 3, "Found: ", c);
                    # Tutaj logika dodawania wyniku
                    Add(cib, res/g);
                fi;

                # 2. Logika Iteratora (Odometer) - zmiana stanu
                pos := rank;
                
                # Szukamy pozycji do inkrementacji (od prawej do lewej)
                # Jeśli cyfra osiągnęła max (g-1), zerujemy ją i przenosimy (carry) dalej
                while pos > 0 and c[pos] = g - 1 do
                    # RESET POZYCJI:
                    # Skoro c[pos] było g-1, to do sumy dodaliśmy (g-1) razy wektor b[pos].
                    # Musimy to teraz odjąć, aby wrócić do stanu, jakby c[pos] było 0.
                    curr_l := curr_l - (g-1) * b[pos];
                    
                    c[pos] := 0;
                    pos := pos - 1;
                od;

                # Jeśli wyszliśmy poza zakres (pos=0), to przeszliśmy cały produkt kartezjański
                if pos = 0 then 
                    break; 
                fi;

                # INKREMENTACJA:
                # Zwiększamy licznik na znalezionej pozycji i aktualizujemy sumę wektorową
                c[pos] := c[pos] + 1;
                curr_l := curr_l + b[pos];
            od;
		od;
		return cib;
    end
);

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
		local res;
		res := [];
		repeat
			if IsGroupByCIBVectorSystemContextRaw(ctx!.data) then
				Add( res, CIBVectorSystemContextRawMat(ctx!.data, 0) / Exponent(ctx) );
			fi;
		until not CIBVectorSystemContextRawInc(ctx!.data);
		return res;
	end
);

CofiniteIntegralBraceVectorSystemsByContext := function(agrp)
	local exp, cib, res, g, name, ctx, rank;

	exp := BraceMaxExponents( agrp );

	cib := [];

	if HasName(agrp) then
		name := Name(agrp);
	else
		name := "anonymous group";
	fi;

	rank := Size(CoboundaryBasisInt(agrp));

	for g in exp do
		Info(InfoCIB, 2, name, ": ", g^rank, " items to process");
		ctx := CIBVectorSystemContext( agrp, g );
		if ctx = fail then
			continue;
		fi;
		res := CofiniteIntegralBraceVectorSystems( ctx );
		Info(InfoCIB, 2, "name: ", name, ", exponent: ", g, " found ", Size(res), " vector systems.");
		Append(cib, res );
	od;
	SetCofiniteIntegralBraceVectorSystems( agrp, cib );
	return cib;
end;

InstallMethod(VectorSystem,
    "for affine crystallographic groups acting on right",
    [IsAffineCrystGroupOnRight and IsStandardAffineCrystGroup],
    function(agrp)
		local pg, ph, d;

		pg := PointGroup(agrp);
		ph := PointHomomorphism(agrp);
		d  := DegreeOfMatrixGroup(pg);

		return List( SSortedList(pg), x->List(PreImagesRepresentative( ph, x )[d+1]{[1..d]}, FractionModOne) );
    end
);

InstallMethod(VectorSystem,
    "for affine crystallographic groups acting on left",
    [IsAffineCrystGroupOnLeft and IsStandardAffineCrystGroup],
	function(agrp)
		return VectorSystem( TransposedMatrixGroup(agrp) );
	end
);
