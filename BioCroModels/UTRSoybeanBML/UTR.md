## Overview

The Thornley utilization-transport-resistance (UTR) model represents the flow of
carbon between plant tissues and its use for growth. Descriptions of several
variations of this model can be found in [Thornley 1972], [Thornley 1998], and
[Thornley 2011].

This repository contains a collection of modules that implement a version of the
UTR model that includes respiration and senescence, but does not (currently)
consider nitrogen. See details below.

## Types of carbon

The following are descriptions of the types of carbon included in the model:

> Here a two-state variable model (structure and substrate; Warren Wilson, 1967)
> is described. 'Substrate' includes storage (easily remobilizable) components.
>
> ([Thornley 2011])

> The experimental scientist will wish to know where this substrate is and how
> it might be measured. In the present model it is envisaged that in each part
> of the plant there is _effectively_ a single pool of substrate which is
> readily available for the purposes of utilization or transport. It is possible
> that a measurement of the sugar concentratration in the sap might provide an
> initially adequate measure of substrate concentration, although in the longer
> term these quantities need to be more precisely defined.
>
> ([Thornley 1972])

> Plant dry matter is considered to consist of structure (X) and of substrates
> (carbon (C) and nitrogen (N)).
>
> ([Thornley 1998])

Note that although the model can include nitrogen, our implementation only
includes carbon. Also note that while Thornley provides several definitions for
substrate carbon, the structural carbon is not usually explicitly defined.
However, it is safe to assume from the definitions above that the structural
carbon must include any components that are not easily remobilizable, such as
carbon incorporated into cell walls or other such components.

Thornley also notes that carbon is not the only component of structural plant
material; for example, he writes that "We assumed that 40% of total dry mass is
structural C" ([Thornley 2011]). Thus, in our model we also distinguish between
structural carbon and biomass.

In this document, we use the following symbols to refer to these different types
of carbon:
- `S` is the mass of substrate carbon
- `X` is the mass of structural carbon
- `M` is the biomass

Subscripts can be used to refer to particular tissues; for example, `S_l` is the
mass of leaf substrate carbon, `X_s` is the mass of stem structural carbon, and
`M_r` is the total root mass.

**Question**: Is biomass proportional to structural carbon or to the sum of
structural and substrate carbon? In other words, is `B ~ X` or `B ~ (X + S)`?

## Processes

The major processes included in Thornley's model are utilization and transport.
However, other processes are also sometimes included, such as senescence (or
litter formation) and respiration. Of course, it is also necessary to consider
photosynthesis as the source of substrate carbon in the leaf. Each of these
processes applies to different types of carbon or to different tissue
components.

> There are only two significant types of process in the plant: transport, and
> chemical/biochemical conversion. Both are necessary and sufficient to
> accomplish allocation. Allocation is the outcome of the processes of substrate
> supply, transport and utilization.
>
> ([Thornley 1998])

> Substrate either remains in the pool of substrate, or it is used for the
> construction of new plant material... The consumption of substrate for growth
> depends upon the local substrate concentration, and obeys a simple relation
> which is widely applied in biochemistry to enzyme-dependent reactions.
>
> ([Thornley 1972])

> Structural dry matter is produced by growth and lost to litter.
>
> ([Thornley 1998])

> It will be assumed that the rate of transport of substrate can be described by
> an expression of the form (`rate of transport of substrate`) =
> (`concentration difference`) / (`a resistance`). Thus, the substrate flows
> down the concentration gradient at a rate proportional to the gradient.
>
> ([Thornley 1972])

To summarize:
- Substrate carbon is subject to transport and utilization.
- Utilization refers to the biochemical conversion of substrate carbon into
  structural carbon.
- Of the utilized substrate carbon, most contributes to the growth of structural
  carbon but some is lost to growth respiration.
- Substrate carbon is transported among plant components, where the transport
  is driven by concentration gradients between the components.
- Structural carbon increases through utilization and decreases through
  senescence.
- Of the structural carbon lost to senescence, most is converted to litter but
  some is recycled into substrate carbon.
- Photosynthesis is a source of substrate carbon and substrate carbon is a sink
  of substrate carbon.

**Note**: Senescence applies to structural carbon but not substrate carbon. If
senescence is also considered to apply to biomass, this would only be consistent
if biomass is proportional to structural carbon alone. (See question above.)

## Equations

### Substrate concentrations

The equations for the utilization and transport rates depend on the
concentrations of substrate carbon in each plant tissue compartment. In
principle, the substrate carbon concentration `s_t` in tissue `t` can be
expressed as `s_t = S_t / V_t`, where `V_t` is the volume of the compartment.
The volume is not generally known, but instead is expected to be proportional to
the mass of structural carbon:

> Assuming that the volume of the reaction in which substrate utilization takes
> place is proportional to the mass of structural C
>
> ([Thornley 2011])

Thus, in our notation, `s_t ~ S_t / X_t`. If densities are similar across plant
tissue compartments, then differences in the "mass ratios" `S / X` should be
similar to differences in the true substrate concentrations.

**Note**: If `B ~ X`, then the ratio `S / B` can also be used. But if the
biomass is not directly proportional to structural carbon, then `S / B` and
`S / X` are not equivalent measures of substrate concentration.

### Utilization

### Transport

### Senescence

### Overall rates of change of substrate and structural carbon

## References

- [Thornley, J. H. M.][Thornley 1972] (1972) "A Model to Describe the
  Partitioning of Photosynthate during Vegetative Plant Growth." _Ann Bot_
  **36**: 419–430.

- [Thornley, J. H. M.][Thornley 1998] (1998) "Modelling Shoot:Root Relations:
  the Only Way Forward?" _Ann Bot_ **81**: 165-171.

- [Thornley, J. H. M][Thornley 2011] (2011) "Plant growth and respiration
  re-visited: maintenance respiration defined – it is an emergent property of,
  not a separate process within, the system – and why the
  respiration: photosynthesis ratio is conservative." _Ann Bot_ **108**:
  1365–1380.

[Thornley 1972]:https://doi.org/10.1093/oxfordjournals.aob.a084601
[Thornley 1998]:https://doi.org/10.1006/anbo.1997.0529
[Thornley 2011]:https://doi.org/10.1093/aob/mcr238
