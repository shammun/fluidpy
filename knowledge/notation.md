# Notation register — symbols used in fluidpy code, in our own words

Filled chapter by chapter by the knowledge-keeper from the symbols the code actually uses (not a copy of the book's
nomenclature list). A symbol whose meaning changes between chapters gets a ⚠️ and one row per meaning.

| Symbol | Meaning | SI unit | Convention / sign | Chapters | Code name |
|---|---|---|---|---|---|
| Γ ⚠️ | Lapse rate: vertical temperature gradient of the environment | K/m (quoted in K/km) | **Ours = Kundu, §1.10: Γ ≡ dT/dz.** Negative when T falls with height. | ch01 → | `dT_dz`, `lapse_rate()` |
| Γ_a ⚠️ | Adiabatic lapse rate: dT/dz of a parcel moved isentropically, Eq. (1.30) Γ_a = −gαT/C_p (= −g/C_p for a perfect gas) | K/m | **Ours (Γ ≡ dT/dz):** Γ_a ≈ −9.8 K/km in dry air; **stable when Γ > Γ_a** (equivalently dθ/dz > 0). | ch01 → | `adiabatic_lapse_rate()` |

**⚠️ Sign trap: two conventions for the lapse rate.** The physics is identical; the sign and the direction of the
stability inequality flip.

| Convention | Definition | Γ_a (dry air) | Stable when | Used by |
|---|---|---|---|---|
| Kundu (ours) | Γ ≡ dT/dz | ≈ −9.8 K/km | Γ > Γ_a | this book, all fluidpy code and explainers |
| Meteorology (standard) | Γ ≡ −dT/dz | ≈ +9.8 K/km | Γ < Γ_a | most atmospheric-science texts, incl. the lapse-rate-feedback literature |

Converting: Γ_met = −Γ_Kundu. When reading a meteorology paper, negate its Γ before using fluidpy, and flip the
inequality. `adiabatic_lapse_rate()` states its convention in the first line of its docstring; the `parcel_stability`
explainer shows the criterion in the form it computes (dT/dz > Γ_a) with the inequality visible.

## Coordinate and sign conventions per chapter
| Chapter | Axes (which is "up") | Origin / reference level | Stress / pressure sign | Reference scales (L, U, T) | Dimensional or non-dimensional code |
|---|---|---|---|---|---|
