def yieldpalsauto():
    n = 25
    species = [i for i in range(n)]
    # first 0 to n_a: spectral
    n_a = 7
    print(f"this DIMet version supports isotopologues "
          f"coloring until m+{n} species")
    spectralPal = sns.color_palette("Spectral", n_a)
    spectralPal = list(spectralPal)[::-1]
    palsautoD = dict()
    k = 0
    while k < n_a:
        palsautoD[species[k]] = spectralPal[k]
        k += 1
    # last until 24
    n_b = n
    addedpal = sns.color_palette("PuOr", 25-n_a)
    addedpal = list(addedpal)
    j = 0
    while k < n_b:  # k starts in 12
        palsautoD[species[k]] = addedpal[j]
        k += 1
        j += 1
    return palsautoD
