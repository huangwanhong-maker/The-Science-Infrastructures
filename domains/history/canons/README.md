# canons: The History Domain's Records

*Status: seed.* A canon is a record of a history's classical facts, in the format of the [codebook](../codebook.md) (§2: JSON lines, one assertion each). A canon is the coarse resolution at which a history is first recorded. Denser sources refine a region of it later without rebuilding the rest. What is built from a canon is the spacetime of *remembered* history: the canon's own selection (the famous, the diplomatic, the written) is visible in its structure, and must be stated with every reading of it.

## 1. The Canons

| Canon | Period | Beings, occasions, ties | Disciplines of its ties | Status | Drafted |
|---|---|---|---|---|---|
| [athens-470-322bc](athens-470-322bc.jsonl) | 470–322 BC | 30, 75, 48 | philosophy 28, relations of states 14, art 2, letters 2, science 2 | seed: not yet checked against its sources | by a language model (Claude), 24 October 2026 |
| [italy-1401-1564](italy-1401-1564.jsonl) | 1401–1564 | 39, 98, 71 | art 48, relations of states 14, philosophy 7, letters 1, science 1 | seed: not yet checked against its sources | by a language model (Claude), 24 October 2026 |
| [science-1543-1727](science-1543-1727.jsonl) | 1543–1727 | 33, 100, 58 | science 53, philosophy 5 | seed: not yet checked against its sources | by a language model (Claude), 24 October 2026 |
| [europe-1785-1828](europe-1785-1828.jsonl) | 1785–1828 | 33, 92, 57 | relations of states 26, music 12, science 7, art 7, letters 3, philosophy 2 | seed: not yet checked against its sources | by a language model (Claude), 24 October 2026 |
| [europe-1871-1914](europe-1871-1914.jsonl) | 1871–1914 | 38, 99, 76 | relations of states 76 | seed: not yet checked against its sources | by a language model (Claude), 24 October 2026 |
| [vienna-1890-1938](vienna-1890-1938.jsonl) | 1890–1938 | 50, 120, 88 | philosophy 23, music 21, art 17, letters 11, science 8, relations of states 8 | seed: not yet checked against its sources | by a language model (Claude), 24 October 2026 |

**Athens, 470–322 BC** runs from the boyhood of Socrates to the deaths of Aristotle and Demosthenes: philosophy foremost (Socrates and his pupils, the sophists, the Academy and the Lyceum), with the politics and wars that framed it, the theatre that mocked it, and the art of Pericles' Athens. The city itself is a being. Its dates are astronomical years of the Julian calendar (399 BC is −0398), most of them a year or two years, as the Athenian year ran from midsummer. Its selection is that of the surviving literature, and some occasions rest on late anecdote, which their confidences say.

**Florence and Rome, 1401–1564** runs from the competition for the Baptistery doors to Michelangelo's death: art foremost, with its patrons, the politics of the Medici, the republic and the popes, and the letters and philosophy around them. Its selection is the canon of art history since Vasari. Its dates are Julian, as the sources write them, and the Florentine year, which began on 25 March, is given in the modern year.

**The Scientific Revolution, 1543–1727** runs from De revolutionibus and the death of Copernicus to the death of Newton: science foremost, with its philosophy, its patrons and its politics; the correspondences of Mersenne and Oldenburg, the trials of Galileo, the founding of the Royal Society and the Académie, and the quarrels over priority. Dates are written as their sources write them: England in Old Style, the Catholic continent in the Gregorian calendar from 1582, Protestant Germany and Denmark in the Julian until 1700.

**Europe, 1785–1828** joins music, letters, philosophy, science, art and the relations of states where they meet. It runs from Mozart's quartets for Haydn to Schubert's death, through the wars of Napoleon and the Congress of Vienna. It was the first canon, drafted to test the method.

**The European alliances, 1871–1914** are the relations of the great powers and a few Balkan states, from the peace of Frankfurt to the declarations of war, with states as collective beings beside their sovereigns and ministers. It is the classic case of structural balance (Antal, Krapivsky and Redner 2005). Its signs are coded occasion by occasion from what the states did, not chosen to make the network balance.

**Vienna, 1890–1938** is the integrated history of one city, from the schoolboy Hofmannsthal reading his verse to Schnitzler to Freud's departure for London. It follows the canonical intellectual history of Vienna 1900 (Schorske; Janik and Toulmin; Kandel), and so the famous, the male, the bourgeois and the written. The Vienna of workers, migrants, most women and Red Vienna's politics is almost absent.

Every canon was drafted from general knowledge. Every assertion names the standard reference it is to be checked against, and none has been checked yet. Each canon's selection is stated in its header. The data library checks a canon's beings against Wikidata (`data/products/history_check.py`). The signs are one coder's.

## 2. Joined Canons

Canons that share beings, occasions or ties are glued along them by `history::join_records` (the theory note, Definition 1.3). A canon that means a being another canon names uses the same id for it. An id identifies across canons only where they agree about what it names: the same Wikidata item, or, where either gives none, the same name. Otherwise the later canon's is renamed `id@canon`, as `alexander` is: Alexander the Great in Athens, Alexander I of Russia in Europe 1785–1828. Of the canons above, the alliances and Vienna share Franz Joseph and Austria-Hungary, and so form one component of the joined history; the others stand apart until a being crosses between them. The History Complex offers the join of every canon as *All canons, joined*.

The facts a canon asserts are not anyone's property; the canon's text is MIT-licensed with the repository.
