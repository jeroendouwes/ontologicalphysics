# Diepgaand Onderzoek: Flyby-anomalie en Kosmische Versnelling

> Onderzoeksdatum: 23 maart 2026
> Status: werkdocument

---

## Deel 1: De Flyby-anomalie — Diepgaand Onderzoek

### 1.1 Huidige status (2020-2026)

De flyby-anomalie is per maart 2026 **nog steeds onverklaard**. Er is geen
wetenschappelijke consensus over de oorzaak. De anomalie werd voor het eerst
opgemerkt bij de Galileo I aardpassage in 1990, en sindsdien bij meerdere
missies waargenomen.

**Compleet overzicht van alle aardpassages:**

| Missie | Datum | Perigeum (km) | v_inf (km/s) | dv (mm/s) | Anomaal? |
|--------|-------|---------------|-------------|-----------|----------|
| Galileo I | dec 1990 | 960 | 8.949 | +3.92 | Ja |
| Galileo II | dec 1992 | 303 | 8.877 | -4.60 | Ja |
| NEAR | jan 1998 | 539 | 6.851 | +13.46 | Ja |
| Cassini | aug 1999 | 1175 | 16.010 | -2 | Ja |
| Rosetta I | mrt 2005 | 1956 | 3.863 | +1.82 | Ja |
| MESSENGER | aug 2005 | 2347 | — | +0.02 | Nee (binnen fout) |
| Rosetta II | nov 2007 | 5322 | — | ~0 | Nee |
| Rosetta III | nov 2009 | 2483 | — | ~0 | Nee |
| Juno | okt 2013 | 559 | — | ~0 | Nee |
| OSIRIS-REx | sep 2017 | — | — | ~0 | Nee |
| BepiColombo | apr 2020 | — | — | ~0 | Nee |

**Recente ontwikkelingen (2024-2026):**

- Een paper uit november 2024 (arXiv:2411.12053) past het variatieprincipe
  van de ART toe met de Lense-Thirring-metriek en claimt de Anderson-formule
  te reproduceren. De auteurs diagonaliseren de Lense-Thirring-metriek om
  de bewegingsvergelijkingen te vinden en voorspellen snelheidsafwijkingen
  van de orde mm/s voor polaire flybys. Dit is nog niet algemeen geaccepteerd.

- In 2019 publiceerde Mirza (MNRAS 489, 3232) een model waarin
  gravitomagnetische frame-dragging wordt versterkt door koppeling met het
  aardmagnetisch veld. Dit zou de anomalie met "ordes van grootte" versterken.
  Echter: Guruprasad (arXiv:1911.05453) bekritiseerde dit en stelde dat
  onafhankelijke radardata met grotere discrepanties, die waren weggelaten
  uit JPL's samenvatting van 2008, een werkelijke bewegingsverandering
  uitsluiten — relativistisch of anderszins.

- Geen enkel voorstel is breed geaccepteerd. Sommige onderzoekers suggereren
  dat verbeterde zwaartekrachtmodellen en nauwkeurigere tracking de schijnbare
  anomalie bij recente missies elimineren, wat erop zou kunnen wijzen dat de
  eerdere waarnemingen systematische fouten bevatten.

**Toekomstige tests:**
- Europa Clipper: aardpassage december 2026
- JUICE (ESA): aardpassage september 2026
- Geen specifieke missie voor de flyby-anomalie gepland of gefinancierd

### 1.2 De Anderson-formule en de aardrotatie

Anderson et al. (2008) vonden empirisch:

    dv/v = K · (cos(delta_in) - cos(delta_out))

met K = 2 · omega_E · R_E / c ≈ 3.099 · 10^-6

waarbij:
- omega_E = 7.292 · 10^-5 rad/s (draaisnelheid Aarde)
- R_E = 6.371 · 10^6 m (straal Aarde)
- delta_in, delta_out = declinatie (breedtegraad) van asymptotische baan
- c = 2.998 · 10^8 m/s

Check: K = 2 · 7.292·10^-5 · 6.371·10^6 / 2.998·10^8 = 3.10 · 10^-6  ✓

De formule past de meeste anomale flybys met ~10% nauwkeurigheid. De
cruciale observatie is dat K de aardrotatie en -straal bevat, en de
anomalie correleert met het verschil in inkomende en uitgaande declinatie.

**Analyse van de structuur:**

De factor omega_E · R_E / c is de ratio van de equatoriale snelheid van
de Aarde (v_eq = 465 m/s) tot de lichtsnelheid:

    v_eq / c = 465 / 3·10^8 = 1.55 · 10^-6

De Anderson-formule zegt dus: dv/v ≈ 2 · (v_eq/c) · (cos(delta_in) - cos(delta_out))

Dit lijkt sterk op een relativistisch effect dat te maken heeft met de
rotatie van het referentiekader. De factor v/c suggereert een eerste-orde
relativistisch effect van de Aardrotatie.

### 1.3 Vergelijking met ART frame-dragging: grootte-orde analyse

**Berekening: Lense-Thirring snelheidsverandering voor een flyby op 500 km hoogte**

Gegeven:
- Aarde impulsmoment: J_E = 7.05 · 10^33 kg·m^2/s
- G = 6.674 · 10^-11 m^3/(kg·s^2)
- c = 2.998 · 10^8 m/s
- r = R_E + 500 km = 6.371 · 10^6 + 5 · 10^5 = 6.871 · 10^6 m
- Flyby-duur typisch ~1000 s (rond perigeum)

De Lense-Thirring precessiesnelheid:

    Omega_LT = 2·G·J / (c^2 · r^3)

Invullen:
    Omega_LT = 2 · 6.674·10^-11 · 7.05·10^33 / ((2.998·10^8)^2 · (6.871·10^6)^3)
             = 9.41·10^23 / (8.988·10^16 · 3.244·10^20)
             = 9.41·10^23 / (2.916·10^37)
             = 3.23 · 10^-14 rad/s

De snelheidsverandering door frame-dragging:

    dv_LT ≈ Omega_LT · r · (typische flyby-duur / iets)

Een meer directe schatting: de frame-dragging snelheid op afstand r is:

    v_FD ≈ Omega_LT · r = 3.23·10^-14 · 6.871·10^6 = 2.22 · 10^-7 m/s
                                                      = 2.22 · 10^-4 mm/s

Dit is de orde van de frame-dragging "windsnelheid" — de snelheid waarmee
het lokale inertieel stelsel meedraait.

Voor een flyby die ~1000 s in het sterke-veld-gebied doorbrengt, is de
cumulatieve snelheidsverandering door Lense-Thirring:

    dv_LT ≈ Omega_LT · r · (effectieve hoek) ~ 10^-7 tot 10^-4 mm/s

**Vergelijking met waarneming:**

| Effect | Grootte-orde |
|--------|-------------|
| Waargenomen anomalie (NEAR) | 13.46 mm/s |
| Waargenomen anomalie (typisch) | 1-10 mm/s |
| ART frame-dragging (Lense-Thirring) | ~10^-4 mm/s |
| Meetonzekerheid | ~0.1 mm/s |

**Conclusie: de standaard ART frame-dragging is ~5 ordes van grootte
te klein om de flyby-anomalie te verklaren.**

De Anderson-formule bevat weliswaar de aardrotatie (omega_E), maar het
effect is ~10^5 keer sterker dan wat Lense-Thirring voorspelt. Als de
anomalie echt is, dan is het NIET standaard frame-dragging.

### 1.4 Kan een ander frame-dragging-model (zoals ORT) een groter effect geven?

**Korte antwoord: nee, niet zonder meer.**

De ORT reproduceert dezelfde Lense-Thirring-formule als de ART in het
zwakke veld (zie frame_dragging_ort.md). De afleiding via ruimtelijke
uitrekking + behoud van impulsmoment geeft:

    Omega_LT = 2GJ / (c^2 r^3)

Dit is identiek aan de ART-voorspelling. De ORT zou alleen een groter
effect kunnen geven als:

1. **Hogere-orde correcties** in het ORT-kader anders zijn dan in de ART.
   In het zwakke veld (aarde) zijn deze correcties echter minuscuul
   (r_s/r ≈ 10^-9 voor de Aarde).

2. **Een ander mechanisme** dan Lense-Thirring verantwoordelijk is.
   De Anderson-formule correleert met de aardrotatie, maar het effect is
   ~10^5 keer sterker dan frame-dragging. Dit suggereert dat als de
   anomalie echt is, er een ANDER mechanisme nodig is.

**Mogelijke ORT-specifieke overwegingen:**

De ORT beschrijft zwaartekracht als snelheidscomponenten (v_grav) in
plaats van metriekkromming. Bij een flyby verandert de richting van
v_grav snel. De combinatie van:
- v_grav als vector (radiaal naar het aardcentrum)
- De draaiing van de Aarde (die de v_grav-richting beïnvloedt via aberratie)
- De asymmetrische geometrie (inkomend vs. uitgaand declinatie)

zou in principe een effect kunnen geven dat verschilt van de standaard
ART-voorspelling. Maar: dit is niet uitgewerkt. Er bestaat geen
ORT-berekening die de Anderson-formule reproduceert of de grootte-orde
van het effect verklaart.

**Eerlijkheid:** De ORT biedt geen verklaring voor de flyby-anomalie.
De standaard frame-dragging is 5 ordes van grootte te klein, en de ORT
reproduceert diezelfde standaard frame-dragging in het zwakke veld.

### 1.5 Waarom vertoonde Juno (2013) geen anomalie?

Juno's aardpassage op 9 oktober 2013 had een perigeum van 559 km — vergelijkbaar
met NEAR (539 km) die de GROOTSTE anomalie vertoonde (+13.46 mm/s). Toch was
er bij Juno geen meetbare anomalie.

**Mogelijke verklaringen:**

1. **Verbeterde modellering:** Tussen 1998 (NEAR) en 2013 (Juno) zijn de
   aardse zwaartekrachtmodellen aanzienlijk verbeterd (GRACE, GOCE). Effecten
   die eerder als "anomaal" werden geïdentificeerd, kunnen bij betere
   modellering verdwijnen. Dit is de meest conservatieve verklaring.

2. **Andere baangeometrie:** De Anderson-formule voorspelt dat de anomalie
   afhangt van cos(delta_in) - cos(delta_out). Als Juno's inkomende en
   uitgaande declinatie vergelijkbaar waren (symmetrische flyby), voorspelt
   de formule dv ≈ 0. De specifieke baangeometrie van Juno kan een
   nul-voorspelling geven.

3. **Missiespecifieke parameters:** De anomalie zou afhankelijk kunnen zijn
   van factoren die per missie variëren: spacecraft-configuratie, antenne-
   geometrie, thermische emissie, atmosferische correcties.

4. **De anomalie is een artefact:** De recente trend (MESSENGER, Rosetta II/III,
   Juno, OSIRIS-REx, BepiColombo — allemaal geen anomalie) suggereert dat
   verbeterde meettechnieken en modellering de schijnbare anomalie elimineren.
   De eerdere anomalieën zouden systematische fouten geweest kunnen zijn.

**De nulresultaten zijn minstens zo informatief als de detecties.** Ze sluiten
elke verklaring uit die een universeel effect voorspelt bij alle aardpassages,
ongeacht baangeometrie.

### 1.6 Verband met de snelheid ten opzichte van de aardrotatie-as?

De Anderson-formule bevat expliciet de declinatie (breedtegraad) van de
asymptotische baan. De declinatie bepaalt de hoek ten opzichte van de
evenaar — en dus ten opzichte van de aardrotatie-as.

**Analyse:**

- cos(delta) = 1 voor delta = 0 (equatoriaal) → maximale bijdrage
- cos(delta) = 0 voor delta = 90° (polair) → geen bijdrage

Een flyby die binnenkomt vanuit de evenaar en vertrekt over de pool heeft
maximale cos(delta_in) - cos(delta_out) = 1 - 0 = 1, en dus maximale anomalie.

Een symmetrische flyby (delta_in ≈ delta_out) geeft cos(delta_in) - cos(delta_out) ≈ 0,
en dus geen anomalie. Dit zou het Juno-nulresultaat kunnen verklaren als de
baangeometrie voldoende symmetrisch was.

**Verband met aardrotatie:** De formule koppelt het effect direct aan de
rotatieparameters van de Aarde (omega_E, R_E). Dit is suggestief voor een
frame-dragging-achtig effect, maar de grootte is ~10^5 keer sterker dan
standaard Lense-Thirring.

De koppeling aan de rotatie-as is een sterke aanwijzing dat de anomalie —
als ze echt is — iets te maken heeft met het roterende referentiekader.
Maar het mechanisme is onbekend.

### 1.7 Samenvatting flyby-anomalie

1. **Status:** Onverklaard per 2026. Recente missies tonen geen anomalie,
   wat de vraag oproept of eerdere detecties systematische fouten waren.
2. **Anderson-formule:** Empirisch, correleert met aardrotatie en declinatie.
   Geen fysieke basis gevonden.
3. **Frame-dragging:** Standaard Lense-Thirring is ~10^5 keer te klein.
   Een versterkt effect (bv. via EM-koppeling, Mirza 2019) is bekritiseerd.
4. **ORT:** Reproduceert dezelfde Lense-Thirring als ART in zwak veld.
   Geen uitgewerkte alternatieve verklaring. Geen groter effect verwacht.
5. **Juno-nulresultaat:** Sluit universele verklaringen uit. Consistent
   met symmetrische baangeometrie of verbeterde modellering.
6. **Toekomst:** Europa Clipper en JUICE flybys in 2026 bieden nieuwe data.

---

## Deel 2: Kosmische Versnelling — Past de ORT bij de waarnemingen?

### 2.1 DESI DR2 resultaten (maart 2025)

Het Dark Energy Spectroscopic Instrument (DESI) publiceerde in maart 2025
de Data Release 2 (DR2), gebaseerd op 14+ miljoen sterrenstelsels en quasars.
Dit zijn de meest precieze baryon-akoestische-oscillatie (BAO) metingen ooit.

**Kernresultaten:**

Het w0waCDM-model (dynamische donkere energie) wordt gefit als:

    w(a) = w0 + wa · (1 - a)

waarbij a = 1/(1+z) de schaalfactor is.

De beste fits geven:
- w0 > -1 (minder negatief dan Lambda)
- wa < 0 (w was negatiever in het verleden)
- w kruist -1 ergens rond z ≈ 0.5

**Statistische significantie:**
- DESI DR2 alleen: LCDM consistent op ~1.5 sigma
- DESI DR2 + CMB: voorkeur voor dynamische w boven LCDM, >2 sigma
- DESI DR2 + CMB + supernovae: 2.8-4.2 sigma (afhankelijk van SN-dataset)

**Kanttekeningen:**
- De significantie hangt sterk af van welke supernovadataset wordt gebruikt
- Er zijn spanningen TUSSEN de BAO-, CMB- en SN-datasets onderling
- Sommige analyses (arXiv:2504.15222, EPJC 2025) betogen dat de claim van
  dynamische donkere energie niet robuust is door deze onderlinge spanningen
- Het is nog NIET op ontdekkingsniveau (5 sigma)

**Wat dit betekent:** Er is een suggestie, maar geen bewijs, dat donkere
energie niet constant is. Als w echt varieert met roodverschuiving, dan
is de kosmologische constante Lambda niet de juiste beschrijving en is er
een nieuw mechanisme nodig.

### 2.2 Wat moet verklaard worden?

De waarneming die "kosmische versnelling" heet, is primair gebaseerd op
de helderheid-afstand relatie van Type Ia supernovae.

**De luminosity-distance relatie:**

In een expanderend heelal geldt:

    d_L(z) = (1+z) · integral_0^z [c · dz' / H(z')]

waarbij H(z) de Hubble-parameter is als functie van roodverschuiving.

Type Ia supernovae zijn "standaardkaarsen" — hun intrinsieke helderheid
is bekend. De waargenomen helderheid geeft de luminosity-distance d_L.
De roodverschuiving z is direct meetbaar.

**De ontdekking (1998):** Verre supernovae (z > 0.5) zijn lichtzwakker
(= verder weg) dan verwacht in een vertragend heelal. Conclusie: de
expansie versnelt.

In LCDM wordt dit verklaard door een kosmologische constante Lambda met
w = -1 (constant). De DESI-hint suggereert w ≠ -1 en variërend.

### 2.3 Hoe beïnvloedt een variërende c_local de observabelen?

De ORT-kosmologie (§12.22-12.25) stelt dat ons heelal de binnenkant van
een zwart gat is, met c_local die varieert als functie van kosmologische
positie en tijd:

    c_local(d, t) = c · sqrt(1 - H(t)^2 · d^2 / c^2)        (formule 78)

**Effect op de luminosity-distance relatie:**

Als c lokaal varieert, verandert de lichtreistijd en dus de afstandsmeting.
De luminosity-distance wordt:

    d_L(z) = (1+z) · integral_0^z [c_local(z') · dz' / H(z')]

In de standaard kosmologie is c_local = c (constant). In de ORT daalt
c_local voor objecten dichtbij de kosmologische horizon (H·d → c).

**Richting van het effect:** Als c_local < c voor verre objecten, dan:

1. Licht reist langzamer → het doet er langer over om ons te bereiken
2. Dit vergroot de effectieve afstand → objecten lijken verder weg
3. Supernovae lijken lichtzwakker → dit IMITEERT versnelde expansie

Dit is precies dezelfde richting als het waargenomen effect. Een dalende
c_local kan in principe versnelde expansie SIMULEREN zonder dat de
expansie daadwerkelijk versnelt.

**Dit is geen nieuw inzicht.** De VSL (Varying Speed of Light) kosmologie
is uitgebreid bestudeerd door onder anderen Albrecht en Magueijo (1999),
Barrow (1999), en recenter door Lee (2021, "meVSL" model). Het kernresultaat:
een dalende c(t) in een Einstein-de Sitter heelal (geen Lambda) kan de
Type Ia supernova-data fitten vergelijkbaar met LCDM.

**Effect op de schijnbare versnelling:**

In VSL-kosmologie wordt de deceleratieparameter q effectief verschoven:

    q_eff = q_intrinsiek - correctieterm(dc/dt)

Als dc/dt < 0 (c daalt), wordt q_eff negatiever (meer versnelling).
Negatieve waarden van de correctie-parameter leiden tot meer schijnbare
versnelling; positieve waarden tot vertraging.

**Effect op het CMB-vermogenspectrum:**

Het CMB-spectrum is extreem gevoelig voor de kosmologische parameters.
Een variërende c(t) beïnvloedt:

1. De geluidshorison bij recombinatie (z ≈ 1100)
2. De Silk-demping-schaal
3. De piekposities en -hoogtes

In de standaard LCDM zijn de CMB-pieken met uitzonderlijke precisie
gefit (~0.5% op de kosmologische parameters). Een VSL-model dat de
supernovae fit, moet OOK het CMB-spectrum fitten. Dit is een STRENGE
eis die de meeste VSL-modellen niet halen.

**Eerlijkheid:** De ORT heeft geen kwantitatief CMB-model. Zonder een
voorspelling van het CMB-spectrum kan de ORT-kosmologie niet volwaardig
concurreren met LCDM.

### 2.4 Kan een variërende c_local Lambda imiteren?

**Ja, in principe.** De wiskundige mogelijkheid is aangetoond in de
VSL-literatuur. Een dalende c(t) kan de luminosity-distance relatie
reproduceren die in LCDM door Lambda wordt verklaard.

**De ORT-specifieke versie:**

In de ORT is c_local niet ad hoc variabel — het volgt uit de
zwaartekrachtspotentiaal:

    c_local = c · sqrt(1 - r_s/r)    (lokaal, §12.1)
    c_local = c · sqrt(1 - H^2·d^2/c^2)  (kosmologisch, §12.24)

De kosmologische versie is de Friedmann-tegenhanger van de lokale formule.
De variatie van c_local is geen vrije parameter — het wordt volledig
bepaald door H(t) en de afstand d.

**Kwantitatieve voorspelling?**

De ORT voorspelt:

    H(t) = c^3 / (2G · M_BH(t))                              (formule 79)
    Lambda_eff ≈ 3c^6 · Omega_Lambda / (4G^2 · M_BH^2)       (formule 80b)

Met M_BH ≈ 9.2 · 10^52 kg geeft dit de juiste ORDE VAN GROOTTE voor
Lambda. Maar dit is geen onafhankelijke voorspelling — de massa M_BH
wordt gekozen om te passen.

**Effectieve w(z):**

De ORT voorspelt kwalitatief dat w varieert met roodverschuiving, omdat
de BH-groeidynamiek niet constant is:

- Als M_BH groeit: H daalt → w_eff > -1 (minder versnelling)
- Als M_BH constant is: H constant → w_eff = -1 (de Sitter)
- Overgangsperiode: w kruist -1

Dit is kwalitatief consistent met de DESI-hint (w kruist -1 rond z ≈ 0.5).
Maar er is GEEN kwantitatieve w(z)-voorspelling uit de ORT.

### 2.5 De richting van het effect: dalend c_local

**Kernvraag:** Als c_local afneemt over kosmische tijd (zoals het
moeder-zwart-gat groeit), geeft dit dan versnelling of vertraging?

**Analyse stap voor stap:**

1. Het moeder-zwart-gat groeit → r_s groeit → de eventhorizon rekt op
2. Binnenin: H(t) = c^3/(2G·M_BH) → als M_BH groeit, daalt H
3. Dalend H → de expansie vertraagt → VERTRAGING, geen versnelling

Maar:

4. Als M_BH stopt met groeien → H wordt constant
5. Constant H → a(t) ∝ exp(Ht) → exponentiële expansie → VERSNELLING
6. De overgang van groei naar stabilisatie geeft de overgang van
   vertraging naar versnelling

**Samenvatting richting:**
- M_BH groeit → H daalt → vertraging (straling- en materietijdperk)
- M_BH stabiliseert → H constant → versnelling (huidig tijdperk)
- M_BH nagenoeg constant → H nagenoeg constant → quasi-de Sitter

De ORT voorspelt dat de versnelde expansie BEGINT wanneer het
moeder-zwart-gat stopt met significant groeien. Dit is het huidig
tijdperk (z < 1), consistent met de waarneming dat de versnelling
begon rond z ≈ 0.7.

### 2.6 Is er een specifieke ORT-voorspelling die verschilt van LCDM?

**Potentiële onderscheidende voorspellingen:**

1. **w(z) varieert:** In LCDM is w = -1 exact. In de ORT varieert w
   omdat M_BH(t) niet precies constant is. De ORT voorspelt kwalitatief
   w > -1 nu, w < -1 in het verleden (quintom B scenario).
   STATUS: consistent met DESI-hint, maar niet kwantitatief.

2. **Geen cosmological coincidence probleem:** In LCDM is het toeval dat
   Omega_Lambda ≈ Omega_m juist nu. In de ORT volgt dit uit de overgang
   van BH-groei naar stabilisatie — een natuurlijk tijdperk.
   STATUS: kwalitatief aantrekkelijk, niet testbaar.

3. **Geen fine-tuning:** De waargenomen Lambda verschilt 10^120 van de
   QFT-verwachting. In de ORT is Lambda_eff bepaald door M_BH, niet
   door vacuümenergie. De vraag verschuift naar "waarom is M_BH deze
   waarde?" — een minder ernstige vraag.
   STATUS: conceptueel beter, niet kwantitatief testbaar.

4. **Verband inflatie-versnelling:** De ORT stelt dat inflatie en de
   huidige versnelling hetzelfde mechanisme zijn (constant M_BH → constant H).
   Als dit klopt, is er een relatie tussen de inflatie-energieschaal en
   de huidige H0. Dit is in principe testbaar.
   STATUS: niet uitgewerkt.

### 2.7 Eerlijke evaluatie

**Wat de ORT WEL biedt:**

- Een conceptueel kader waarin versnelde expansie VOLGT uit de
  BH-groeidynamiek, zonder Lambda als vrije parameter
- De juiste orde van grootte voor Lambda_eff
- Kwalitatieve consistentie met de DESI-hint van variërend w
- Oplossing van het fine-tuning-probleem (10^120 ordes verdwijnen)
- Unificatie van inflatie en huidige versnelling

**Wat de ORT NIET heeft:**

- **Geen kwantitatief kosmologisch model.** Er is geen ORT-equivalent
  van de Friedmann-vergelijkingen met w(z) als output.
- **Geen CMB-voorspelling.** Zonder een voorspelling van het
  CMB-vermogenspectrum is de ORT-kosmologie niet competitief met LCDM.
- **Geen specifieke d_L(z) voorspelling.** De luminosity-distance relatie
  als functie van roodverschuiving is niet berekend in het ORT-kader.
- **Geen BAO-voorspelling.** De BAO-schaal is niet berekend.
- **Geen structuurvorming.** Hoe vormen sterrenstelsels in de ORT-kosmologie?
- **M_BH is een vrije parameter.** De ORT vervangt Lambda door M_BH,
  maar voorspelt M_BH niet.

**De kern:** De ORT-kosmologie is op dit moment een **kwalitatief raamwerk**
— conceptueel aantrekkelijk, maar zonder kwantitatieve voorspellingen.
Het is geen alternatief voor LCDM in de zin dat het data kan fitten.
Het is een hypothese die richting geeft voor verder onderzoek.

### 2.8 Wat zou nodig zijn?

Om de ORT-kosmologie van hypothese naar testbaar model te brengen:

1. **Afleiden van w(z) uit M_BH(t):** Specificeer hoe M_BH groeit als
   functie van kosmische tijd. De BH-accretiehistorie bepaalt H(t) en
   daarmee w(z). Dit is de meest haalbare stap.

2. **Berekenen van d_L(z):** Met H(t) uit stap 1 kan de luminosity-distance
   relatie berekend worden en vergeleken met supernova-data.

3. **CMB-spectrum berekenen:** Dit vereist de volledige kosmologische
   perturbatietheorie in het ORT-kader — een enorme inspanning.

4. **BAO-voorspelling:** De BAO-schaal hangt af van de geluidshorison
   bij recombinatie, die beïnvloed wordt door c_local(t).

5. **Vergelijking met DESI DR2:** Pas het ORT-model met M_BH(t) aan op
   de DESI BAO-data en vergelijk de chi-kwadraat met LCDM en w0waCDM.

---

## Deel 3: Synthese — Verbanden en Open Vragen

### 3.1 Flyby-anomalie en ORT

De flyby-anomalie is NIET het sterkste testgeval voor de ORT:
- De ORT reproduceert standaard Lense-Thirring (5 ordes te klein)
- De anomalie is mogelijk een meetartefact (recente nulresultaten)
- Er is geen ORT-mechanisme dat het grootte-orde-probleem oplost

De flyby-anomalie is wel een interessante puzzel die met nieuwe data
(Europa Clipper, JUICE 2026) verder onderzocht zal worden.

### 3.2 Kosmische versnelling en ORT

De kosmische versnelling is WEL een relevant testgeval:
- De ORT biedt een conceptueel alternatief voor Lambda
- De richting van het effect (c_local-dynamiek) is correct
- De DESI-hint naar variërend w past kwalitatief
- Maar: geen kwantitatieve voorspellingen

### 3.3 Prioriteiten voor de ORT

1. **Hoogste prioriteit:** Een kwantitatieve w(z) uit M_BH(t). Dit is
   de meest directe test tegen DESI-data.

2. **Tweede prioriteit:** De flyby-anomalie berekenen met het volledige
   ORT-kader (v_grav + frame-dragging + c_local-variatie). Zelfs als het
   antwoord "nul effect" is, is dat informatief.

3. **Lange termijn:** CMB-perturbatietheorie in ORT-kader.

---

## Bronnen

### Flyby-anomalie
- [Flyby anomaly - Grokipedia](https://grokipedia.com/page/Flyby_anomaly)
- [Flyby Anomaly in the Variation Principle of GR (2024)](https://arxiv.org/abs/2411.12053)
- [Flyby anomaly and gravitational-magnetic field frame-dragging (Mirza 2019)](https://academic.oup.com/mnras/article/489/3/3232/5555577)
- [Comment on Mirza 2019 (Guruprasad)](https://arxiv.org/abs/1911.05453)
- [The flyby anomaly: A case for strong gravitomagnetism? (2015)](https://arxiv.org/abs/1505.06884)
- [The flyby anomaly: A multivariate analysis approach (2017)](https://arxiv.org/abs/1701.05735)
- [Flyby anomaly - IFLScience](https://www.iflscience.com/flyby-anomaly-the-unexplained-phenomenon-affecting-several-nasa-spacecraft-76014)

### DESI DR2 / Kosmische versnelling
- [DESI DR2 Results Guide (maart 2025)](https://www.desi.lbl.gov/2025/03/19/desi-dr2-results-march-19-guide/)
- [Key Results from DESI DR2 - Astrobites](https://astrobites.org/2025/10/06/desi-dr2-part1/)
- [Dynamical dark energy in light of DESI DR2 - Nature Astronomy](https://www.nature.com/articles/s41550-025-02669-6)
- [Did DESI DR2 truly reveal dynamical dark energy? - EPJC](https://link.springer.com/article/10.1140/epjc/s10052-025-15076-y)
- [DESI DR2 Results II - arXiv](https://arxiv.org/abs/2503.14738)
- [DESI hints at evolving dark energy - CERN Courier](https://cerncourier.com/desi-hints-at-evolving-dark-energy/)

### VSL (Varying Speed of Light) kosmologie
- [A statefinder luminosity distance formula in VSL cosmology](https://www.researchgate.net/publication/262806045)
- [The minimally extended VSL (meVSL)](https://arxiv.org/abs/2011.09274)
- [A mechanism to generate VSL via Higgs-dilaton coupling (2025)](https://link.springer.com/article/10.1140/epjc/s10052-025-14082-4)
- [Redshift drift and luminosity distance in VSL cosmology](https://www.researchgate.net/publication/258105992)
- [Current constraints on meVSL through cosmic distance duality](https://arxiv.org/html/2505.15768)

### Lense-Thirring / Frame-dragging
- [Lense-Thirring precession - Wikipedia](https://en.wikipedia.org/wiki/Lense%E2%80%93Thirring_precession)
- [Frame-dragging - Wikipedia](https://en.wikipedia.org/wiki/Frame-dragging)
