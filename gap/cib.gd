#
# cib: Cofinite Integral Braces in GAP
#
#! @Chapter Introduction
#!
#! cib is a package which does some
#! interesting and cool things
#!
#! @Chapter Functionality
#!
#!
#! @Section Additional methods for affine crystallographic groups
#!
#! This section will describe the properties and methods, which
#! can be thought of as an extension to the <Package>Cryst</Package> package.

#! @Arguments grp
#! @Returns true if <A>grp</A> is an integral affine crystallographic group.
#! @Description
#! An affine cryst group is called **standard**, if its <Ref BookName="Cryst" Func="InternalBasis"/> forms an identity matrix.
#! An affine cryst group is called **integral**, if it is standard, and its <Ref BookName="Cryst" Func="PointGroup"/> is an integer
#! matrix group.
DeclareProperty("IsIntegralAffineCrystGroup", IsStandardAffineCrystGroup);

#! @Arguments grp
#! @Returns a basis of rational space of vector systems for <A>grp</A>, which define equivalence of extensions with <A>grp</A>.
#! @Description
#! An argument <A>grp</A> is an integral affine crystallographic group.
#! Let $G$ be its <Ref BookName="Cryst" Func="PointGroup"/>. Then $G$ acts naturally on $\mathbb{Z}^n$.
#! With the assumption that <A>grp</A> acts from the right, we store the sorted list of $G$.
#! If <A>grp</A> acts from the left, we store the sorted list of $G^t$, the transpose of $G$.
#! Assume that this sorted list is of the form $(g_1, g_2, \ldots, g_k)$.
#! Then every element of the basis is a $k \times n$ matrix. Let 
#! $$C = ( c_1, c_2, \ldots, c_k )$$ 
#! be any rational linear combination of the basis elements. Then the group of matrices of the form
#! <Alt Only="Text"><Verb>
#!         [ g_i   0 ]
#!         [ c_i+z 1 ],</Verb>
#! </Alt>
#! <Alt Not="Text"><Display>
#! \left[\begin{array}{rr} g_i &amp; 0 \\ c_i +z &amp; 1 \end{array}\right],
#! </Display></Alt>
#! where $1 \leq i \leq k$ and $z \in \mathbb{Z}^n$, is a split extension of $\mathbb{Z}^n$ by $G$.
#! Conversely, every split extension of $\mathbb{Z}^n$ by $G$ is obtained in this way.
#! 
#! **Important note:** the basis is always calculated for the right action, hence if <A>grp</A> acts from the left,
#! then the transposed matrices should be used.
#! <Label Name="CoboundaryBasisInt"/>
DeclareAttribute("CoboundaryBasisInt", IsIntegralAffineCrystGroup);

#! @Arguments grp
#! @Returns a vector system for <A>grp</A>.
#! @Description
#! As for <Ref Attr="CoboundaryBasisInt"/>, the returned list corresponds to the list of sorted
#! elements of the <Ref BookName="Cryst" Func="PointGroup"/> of <A>grp</A>.
#! If <A>grp</A> acts from the left, then the transposed matrices should be used.
DeclareAttribute("VectorSystem", IsStandardAffineCrystGroup);

#! @Section Braces

#! @Arguments grp
#! @Returns the maximal exponents for the additive group of a brace with the multiplicative group isomorphic to <A>grp</A>.
#! @Description
#! In an obvious way, the exponents are bounded by the order of the group <A>grp</A>.
#! However, in many cases the maximal exponents are much smaller.
#! This attribute stores the list of maximal exponents, since there may be more than one divisor of the order of <A>grp</A> involved.
DeclareAttribute("BraceMaxExponents", IsGroup);

#! @Description
#! This function clears the cache of maximal exponents for braces with given abelian groups,
#! and the data that is needed for their computation.
#! Having the cache is useful, since the calculation involves computation of automorphism groups of possibly many abelian groups.
DeclareGlobalFunction("ClearCIBMaxExponentCache");

BindGlobal("CIBMaxExponentCache", rec( abelian_groups := [], exponents := [] ) );

#! @Section Cofinite integral braces

DeclareCategory("IsCIBVectorSystemContext", IsComponentObjectRep and IsAttributeStoringRep);
DeclareCategory("IsCIBVectorSystemContextData", IsObject);

DeclareGlobalFunction("CIBVectorSystemContext");

#! @Arguments obj
#! @Returns all possible vector systems of an affine cryst group attached to <A>obj</A>.
#! @Description
#! As in <Ref Attr="CoboundaryBasisInt"/>, the returned list corresponds to the list of sorted
#! elements of the <Ref BookName="Cryst" Func="PointGroup"/> of the affine crystallographic group attached to <A>obj</A>.
DeclareAttribute("CofiniteIntegralBraceVectorSystems", IsObject);
