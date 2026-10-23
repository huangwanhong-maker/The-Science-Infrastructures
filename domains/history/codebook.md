# The Codebook of the History Domain

*Version 1, 2026-10-24. Status: draft.* The conventions by which a history is recorded, so that two coders who read the same sources write the same record, and a record can be checked against them. A canon names the version of the codebook it follows; a change of the codebook is a change of every record coded by it ([ADR 0041](../../development/decisions/0041-the-history-domain.md)).

## 1. What Is Recorded

A record holds **assertions**, never facts: each says what a source says, with a confidence, and with who drafted it and who checked it. Two sources that disagree are two assertions, and the disagreement is kept (the factoid model of prosopography; Bradley and Short 2005). From a record the domain's construction derives a 2-complex without further judgement (the theory note, [§2](theory/the-history-complex.md)):

- an **occasion** is a vertex;
- a **being** is a worldline: its occasions in order, consecutive ones joined by an edge;
- a **tie** between two beings is a chain of 2-cells, the **lenses**, one between each two consecutive meetings of the two within the tie's life.

Only three kinds of thing are therefore coded: occasions with the beings present, and ties with their signs where they are attested. Spins, amplitudes and every other label of the model are never coded.

## 2. The Format

A record is a file of JSON lines (`.jsonl`, UTF-8), one assertion to a line, in any order. Every line has a `kind`. Any line may carry `note` (prose), `drafted_by` and `checked_by` (a person's name, or `model:<model>` for a language model), and `confidence` where §7 asks for one.

| `kind` | Fields | |
|---|---|---|
| `canon` | `id`, `title`, `about`, `codebook` (this version), `period` (EDTF interval), `status`, `drafted_by`, `drafted` (date) | exactly one, the first line |
| `source` | `id`, `cite` | a reference that assertions name |
| `being` | `id`, `name`, `type` (§3), `lived` (EDTF interval of the years of birth and death), `wikidata` (its QID, when known) | |
| `occasion` | `id`, `date` (§4), `calendar` (§4), `place`, `type` (§5), `label`, `present` (a list of `[being, role]`), `sources` (a list of source ids), `confidence` | a vertex |
| `tie` | `id`, `a`, `b` (beings), `type` (§6), `colour` (§6), `attested` (a list of `[occasion, sign]`), `ends` (an occasion, or null), `sources`, `confidence` | a chain of lenses |

Identifiers are lowercase ASCII words joined by `_`; occasions begin with `o` and the year (`o1808_erfurt`), ties are named by their two beings (`napoleon~alexander`).

## 3. Beings

A being is anything that persists and attends occasions: a **person**, a **group** (a school, a circle, a family), an **institution** (a society, an academy, a court, a journal), a **state**, or a **work** (a book, a symphony, a painting, a treaty's text) when the work itself is carried from occasion to occasion. A relation among more than two beings is held through a being of its own, a school or a work or an alliance, with a tie to each member: the duality of persons and groups (Breiger 1974).

A being's identity is fixed by its Wikidata identifier where there is one, checked against Wikidata by the data library (`data/products/history_check.py`); `lived` serves the check and the eye.

## 4. Dates and Calendars

- **Dates** are EDTF (ISO 8601-2:2019, level 1), in the subset: a year `1795`, a month `1797-03`, a day `1785-02-12`; an interval of two of these `1797-03/1797-05`; a qualifier on the whole date: `?` uncertain, `~` approximate, `%` both (`1795?`, `1790-12~`). A date denotes every day from its first to its last.
- **Calendars.** A date is written as the source writes it, and `calendar` says in which calendar:

  | `calendar` | Meaning |
  |---|---|
  | `gregorian` (the default) | the Gregorian calendar, proleptic before 1582 |
  | `julian` | the Julian calendar, the year beginning on 1 January |
  | `julian-march` | Old Style as England wrote it until 1752: Julian, the legal year beginning on 25 March, so that 1 January to 24 March belong to the year before |

  England and its colonies used Old Style until 1752; the Protestant German states the Julian calendar until 1700; Russia until 1918 (the October Revolution, 25 October 1917 Julian, was 7 November Gregorian). The domain converts every date to Julian day numbers, and orders occasions by them.
- **Order.** A being's occasions are ordered by their first day, then their last, then the order of the file. When two of a being's occasions overlap in time (a congress of weeks and an audience during it), the order rests on the file, and the construction says so: a finer resolution is wanted there.
- **A reader's dates.** The dates a reader gives, not the record (the cut of a moment, the end of a year, a year on the axis of time), name no calendar. They are read in the historians' convention (Meeus 1998, ch. 7): the Julian calendar before the Gregorian reform of 15 October 1582, the Gregorian from then on. A record of antiquity or of the Renaissance, dated in the Julian calendar, then has its years end where the reader's do. Read in the proleptic Gregorian, an occasion of Athens dated only "337 BC" would fall five days before the end of 338 BC.

## 5. Occasions

**The test of an occasion is a date and a list of participants, attested by a source.** A period, a movement, a reign or a war of years is not an occasion: it is a type of ties, a sub-complex, or an occasion of a coarser resolution than the record's.

**Types:** `meeting`, `letter`, `lesson`, `publication`, `dedication`, `performance`, `exhibition`, `commission`, `appointment`, `treaty`, `congress`, `battle`, `ceremony`, `journey`, `death`, `other`.

**Who is present.** Every being whose ties the occasion changes, in its role: the writer and the recipient of a letter; the dedicator and the dedicatee; the reviewer and the reviewed; the signatories and the parties a treaty binds (role `party`); the person a work is about (role `subject`); the target of a secret alliance (role `target`); the dead (`died`) and, derived by the construction, the bereaved (§8). Roles are free words, but these are used where they fit: `host`, `guest`, `author`, `recipient`, `teacher`, `pupil`, `patron`, `composer`, `painter`, `performer`, `audience`, `sovereign`, `minister`, `signatory`, `party`, `commander`, `subject`, `target`, `witness`, `died`, `departing`.

## 6. Ties, Their Colours and Their Signs

A tie is dyadic, between two beings, of a type, and has one **colour**, the discipline of its theme:

| Colour | Its ties |
|---|---|
| `science` | collaboration, correspondence and controversy in natural knowledge; the patronage of science |
| `art` | the visual arts: teaching, commission, patronage, rivalry |
| `music` | teaching, dedication, performance, patronage, rivalry in music |
| `letters` | literature: friendship, collaboration and quarrel among writers; the patronage of letters |
| `philosophy` | teaching, influence, dispute and patronage in philosophy |
| `international-relations` | states, sovereigns, ministers and commanders: alliance, enmity, service, diplomacy |

A tie between beings of two disciplines takes the colour of its theme: Beethoven's dedication of a symphony to Bonaparte is `music`; Goethe's audience with Napoleon is `letters`.

**The sign of a tie at an occasion** records the relation *as it stands at that occasion*, read from what the two do there, not from their motives:

| Type | +1 | −1 |
|---|---|---|
| `alliance` (states, sovereigns, commanders) | on the same side: allied, a peace that joins them | at war; a peace dictated by a victor to the vanquished; one moving against the other |
| `diplomacy` (ministers, envoys) | acting together at the occasion | opposed at the occasion |
| `service` (minister and sovereign) | serving | acting against: betrayal, open opposition, dismissal |
| `teaching` | lessons, dedication, praise | rupture, an offence taken, public disparagement |
| `patronage`, `commission` | support given and taken | refusal, a break |
| `friendship`, `collaboration` | friendship, work together | a quarrel, a break |
| `admiration` (one-sided) | admiration expressed | admiration withdrawn |
| `rivalry` | regard shown | rivalry acted out |
| `kinship`, `marriage` | the bond in good standing | estrangement |

- **Secret acts** count at the occasion where they are done, though the other party does not know them (Talleyrand's counsel to the tsar at Erfurt): the record's frame is the historian's, not a participant's.
- **Mixed acts.** When a source attests acts of both signs at one occasion, the sign is that of the act the occasion is remembered for, and the confidence is lowered.
- **Attest only what is attested.** A tie's sign is written only at occasions where a source says something about the relation; at other meetings of the two it is left to the model (§8).

## 7. Confidence

| Confidence | Meaning |
|---|---|
| 1.0 | documented by a primary source, uncontested (a death, a treaty's signature) |
| 0.9 | the standard references agree |
| 0.8 | standard, with minor uncertainty of date or detail |
| 0.7 | plausible, from one secondary source |
| 0.6 | uncertain: to check before the record is used for anything but a draft |

Nothing below 0.6 enters a canon.

## 8. The Life of a Tie, and the Bereaved

- **A tie lives from its first attestation to its last**, or to its `ends` occasion when one is given: a death, a break, a departure. It never lives beyond what is recorded; a tie the sources stop speaking of ends there, and the record says nothing of it after.
- **The bereaved.** When a tie ends at the death or the departure of one of its beings, and the other is not present at that occasion, the construction adds the other as present in the role `bereaved` (derived, marked so), because the occasion changes its ties.
- **Between attestations** a tie passes through the meetings of its two beings at which nothing is said of it; its sign there is latent, and the readings carry it forward with the persistence of its colour, never as if it were attested.

## 9. Provenance and Checking

- `drafted_by` names who wrote an assertion; a canon drafted by a language model says so on its first line, and its assertions are proposals until checked.
- `checked_by` names the person who checked an assertion against its source.
- The data library checks a record's beings against Wikidata (identity, years of birth and death) and reports agreement, disagreement and what it could not match; a disagreement is resolved by a person, in the record, never by the check.
- **Two coders.** Before a record's signs are used for more than exploration, a second coder codes a sample of its ties by this codebook, and the agreement (Krippendorff's α) is reported with the record.

## 10. References

Bradley, J. and Short, H. (2005), Texts into databases: the evolving field of new-style prosopography, *Literary and Linguistic Computing* 20 suppl., 3–24. Breiger, R. L. (1974), The duality of persons and groups, *Social Forces* 53, 181–190. ISO 8601-2:2019 (EDTF). Krippendorff, K. (2018), *Content Analysis*, 4th ed. Meeus, J. (1998), *Astronomical Algorithms*, 2nd ed., ch. 7, Julian Day. *Explanatory Supplement to the Astronomical Almanac* (3rd ed., 2013), ch. 15, Calendars.
