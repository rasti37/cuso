def print_results(
    problem,
    bounds_guess,
    monomial_vertices,
    precompute_multiplicity,
    tshifts_to_include,
    delta_star,
    tau,
    coprime_condition,
):
    print("We are examining the Coppersmith problem defined by")
    print()
    print(problem)
    print()
    print(f"We run Algorithm 3 of eprint:2024/1577 on this problem with parameters")
    print(f"\tprecomputation multiplicity (k_pre): {precompute_multiplicity}")
    print(f"\tMonomial vertices (M_vert): {monomial_vertices}")
    print(f"\tGuessed bounds (X_guess): {bounds_guess}")
    print(f"\tt-shift mask: {tshifts_to_include}")
    print()
    if len(coprime_condition) > 0:
        print("With these parameters, we can compute the precomputation set so long as "
              f"{problem.modulus_multiple_name} is coprime with")
        for cf in coprime_condition:
            print(f"\t{cf}")
    else:
        print("With these parameters, we can always compute the precomputation set.")
    print()
    print("According to Algorithm 3, this problem is heuristically solvable for delta < delta^* when")
    print(f"\ttau = {tau}")
    print("and")
    print(f"\tdelta^* = {delta_star}")
    print()
    print("Substituting this back into the symbolic bounds, this is")
    p = problem.modulus_name
    N = problem.modulus_multiple_name
    for xi, bound_i in zip(problem.ring.gens(), problem.bounds):
        bound_val = bound_i.constant_term + bound_i.linear_term * delta_star
        s = f"\t{xi} < {p}^" + r"{" + f"{bound_val}" + r"}"
        if problem.divisor_bound != 1:
            s += f" = {N}^" + r"{" + f"{bound_val*problem.divisor_bound}" + r"}"
        print(s)
    print()
