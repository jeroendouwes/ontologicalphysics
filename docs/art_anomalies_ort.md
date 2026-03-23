# Anomalieen en Spanningen in de ART — Analyse vanuit de ORT

> Onderzoek: maart 2026
> Status: werkdocument

---

## Inhoudsopgave

1. [Flyby-anomalie](#1-flyby-anomalie)
2. [Rotatiecurves van melkwegstelsels](#2-rotatiecurves-van-melkwegstelsels)
3. [Hubble-spanning](#3-hubble-spanning)
4. [Kosmische versnelling / donkere energie](#4-kosmische-versnelling--donkere-energie)
5. [Zwaartekrachtsgolf-echo's](#5-zwaartekrachtsgolf-echos)
6. [Pioneer-anomalie](#6-pioneer-anomalie)
7. [Samenvatting en conclusies](#7-samenvatting-en-conclusies)

---

## 1. Flyby-anomalie

### 1.1 De waarneming

Bij meerdere ruimtevluchten die een zwaartekracht-slingermanoeuvre (gravity assist)
rond de Aarde uitvoerden, is een klein maar meetbaar verschil gevonden tussen de
voorspelde en de waargenomen snelheid na het passeren:

| Missie | Datum | Perigeum (km) | Anomalie (mm/s) |
|--------|-------|---------------|-----------------|
| Galileo I | dec 1990 | 960 | +3.92 |
| NEAR | jan 1998 | 539 | +13.46 |
| Cassini | aug 1999 | 1175 | -2 |
| Rosetta I | mrt 2005 | 1956 | +1.82 |
| Messenger | aug 2005 | 2347 | +0.02 |
| Rosetta II | nov 2007 | 5322 | ~0 |
| Juno | okt 2013 | 559 | ~0 (niet detecteerbaar) |

De afwijkingen zijn klein (orde mm/s op snelheden van km/s), maar significant
groter dan de meetfout. Anderson et al. (2008) vonden een empirische formule
die de anomalie correleert met het verschil in inkomende en uitgaande
declinatiehoek:

    dv/v = 2 · omega_E · R_E / c · (cos(delta_in) - cos(delta_out))

waarbij omega_E de draaihoeksnelheid van de Aarde is en delta de declinatie
(breedtegraad) van de asymptotische baan. De anomalie lijkt te correleren met
de Aardrotatie en met asymmetrie in de breedtegraad van de baan.

### 1.2 Hoe de ART het verklaart (of niet)

De ART verklaart dit effect **niet**. De standaard gravitationele effecten
(Schwarzschild + frame-dragging door de Aardrotatie) zijn veel te klein om
de waargenomen afwijkingen te verklaren. De Lense-Thirring-precessie door de
Aardrotatie levert effecten in de orde van nano-m/s, niet milli-m/s.

Systematische fouten (atmosferische weerstand, thermische effecten, Aardgetijden)
zijn uitgebreid onderzocht maar verklaren het patroon niet. Het feit dat sommige
flybys (Juno, Rosetta II) géén anomalie vertonen, maakt het patroon complexer.

**Status (2026):** De flyby-anomalie blijft onverklaard. Een paper uit 2024
(arXiv:2411.12053) claimt dat een variatieprincipe-uitbreiding van de ART in
roterende referentiekaders de Anderson-formule reproduceert, maar dit is nog
niet algemeen geaccepteerd.

### 1.3 Hoe de ORT het mogelijk anders verklaart

De ORT biedt een **potentieel relevant kader**, maar geen uitgewerkte verklaring:

**Wat de ORT wél biedt:**
- v_grav = sqrt(2GM/r) als een concreet snelheidseffect met een richting
  (radiaal, naar de massa toe). Dit is equivalent aan Gullstrand-Painlevé
  coördinaten.
- Frame-dragging in de ORT (§12.27) via het c_local(r,theta)-veld en het
  sleepveld omega(r,theta). De Kerr-uitbreiding geeft theta-afhankelijke
  effecten.
- De isotrope c_local-bijdrage van de Aarde als schil: binnen de Aarde
  verlaagt de schilmassa c_local isotroop (schiltheorema).

**Mogelijke relevantie:**
De Anderson-formule correleert met de Aardrotatie. In de ORT is frame-dragging
een gevolg van ruimtelijke uitrekking + behoud van impulsmoment. Dit geeft
dezelfde Lense-Thirring-formules als de ART in het verre veld, maar zou in
het nabije veld (laag perigeum) subtiel kunnen afwijken doordat de ORT de
effecten beschrijft als snelheidscomponenten in plaats van metriekkromming.

De vraag is of de **combinatie** van:
1. Het radiaal gerichte v_grav-veld
2. Het rotatiegebonden sleepveld omega(r,theta)
3. De isotrope c_local-verlaging door de Aardmassa

een netto-effect geeft dat verschilt van de ART bij asymmetrische flybys.
De grootte-orde van het verschil zou mm/s moeten zijn.

**Eerlijkheid:** Er is geen uitgewerkte ORT-berekening die de flyby-anomalie
reproduceert. De ORT reproduceert de ART-voorspellingen in het zwakke veld,
dus een verschil zou moeten komen uit hogere-orde correcties. Dit is een
**open vraag** en een potentieel interessant testgeval.

### 1.4 Meetprecisie

- Doppler-tracking: ~0.1 mm/s precisie
- De anomalie (1-13 mm/s) is 10-100x groter dan de meetfout
- Het niet-detecteren bij Juno en Rosetta II is even informatief als de detecties

---

## 2. Rotatiecurves van melkwegstelsels

### 2.1 De waarneming

De rotatiesnelheden van sterren en gas in spiraalstelstels blijven constant
(of stijgen zelfs licht) op grote afstanden van het centrum, in plaats van
af te nemen als 1/sqrt(r) zoals verwacht bij een geconcentreerde massa:

    v_Kepler(r) = sqrt(GM(<r) / r)  -->  v zou moeten dalen als 1/sqrt(r)
    v_waargenomen(r)                 -->  blijft constant (~200-300 km/s)

Dit is een van de sterkste aanwijzingen voor het bestaan van donkere materie:
een onzichtbare massacomponent die een halo vormt rond elk sterrenstelsel en
waarvan M(r) proportioneel met r toeneemt op grote afstanden.

### 2.2 Hoe de ART het verklaart

De ART zelf voorspelt dit effect niet. Het wordt verklaard door het toevoegen
van **donkere materie** — een extra component die:
- ~85% van alle materie in het heelal uitmaakt
- niet elektromagnetisch wisselwerkt (onzichtbaar)
- een halo vormt rond sterrenstelsels met dichtheidsprofiel rho proportioneel met 1/r^2
- consistent is met de kosmische achtergrondstraling (CMB), lenzing, en
  structuurvorming op grote schaal

**Alternatieven:**
- MOND (Modified Newtonian Dynamics, Milgrom 1983): past de zwaartekrachtwet
  aan onder een kritische versnelling a_0 ~ 1.2 · 10^-10 m/s^2. Past goed op
  individuele melkwegstelsels maar heeft moeite met clusters en de CMB.
- Recente (2025) Gaia-data toont een dalende rotatiecurve voor de Melkweg,
  wat problematisch is voor MOND.
- Een 2026 paper stelt frame-afhankelijke effecten voor geassocieerd met de
  Melkweg als verklaring (single parameter model).

**Status (2026):** Donkere materie blijft de consensus-verklaring. Directe
detectie-experimenten (LZ, XENONnT, PandaX) hebben nog geen WIMP-deeltjes
gevonden. De bovengrens op de werkzame doorsnede wordt steeds strenger.

### 2.3 Hoe de ORT het mogelijk anders verklaart

De ORT biedt twee potentiele mechanismen (§12.25):

**1. Donkere materie als randeffect (speculatief)**

Als ons heelal de binnenkant van een zwart gat is (§12.22-12.24), dan
ervaren objecten extra zwaartekracht van de "schil" — de rand van het heelal:

    g_eff(r) = g_Newton(r) + g_boundary(r)
    g_boundary ~ c^2/R_H · (r/R_H)

Dit is een lineair groeiende correctie die verwaarloosbaar is op kleine schalen
maar significant wordt op galactische schalen. Kwalitatief is dit precies wat
nodig is: een extra kracht die toeneemt met afstand tot het centrum.

**2. Multi-mass c_local-effect (schiltheorema)**

De ORT maakt een fundamenteel onderscheid:
- **v_grav** (scalair, potentiaal): bepaalt c_local, geen kracht
- **v_netto** (vector, veld): geeft de gravitationele versnelling

Binnen een bolvormige schil geldt: g = 0 (geen nettokracht), maar
c_local < c (lagere potentiaal). Een sterrenstelsel is geen puntvormige massa
maar een uitgebreide verdeling. De c_local op een punt wordt bepaald door de
potentiaal van ALLE omringende massa, niet alleen de massa binnen de baan.

De vraag is: leidt de isotrope c_local-verlaging door de massa buiten de baan
tot een meetbaar effect op de rotatiesnelheid? In de standaard ART is dit
effect nul (de Birkhoff-stelling garandeert dat de buitenschil geen invloed
heeft op de interne dynamica). Maar de Birkhoff-stelling gaat over krachten
(gradienten), niet over potentialen. De ORT beschrijft zwaartekracht als een
potentiaal-effect (c_local), niet primair als een kracht.

**Eerlijkheid:**
- De randcorrectie-formule (81b) is **kwalitatief**, geen kwantitatieve
  voorspelling. Het is niet aangetoond dat dit de juiste rotatieprofiel-vorm
  geeft (constant v(r) vereist een specifieke r-afhankelijkheid).
- De c_local-potentiaal-bijdrage van de buitenschil geeft tijddilatatie maar
  geen kracht (schiltheorema geldt ook in de ORT).
- Er is **geen uitgewerkte ORT-berekening** die specifieke rotatiecurves fit.
- Het MODEL.md zegt dit expliciet: "Dit geeft een kwalitatieve verklaring,
  geen kwantitatieve voorspelling."

### 2.4 Meetprecisie

- HI 21cm rotatiemetingen: ~5-10 km/s precisie
- Optische spectroscopie: ~1-5 km/s
- De afwijking is ~100-200 km/s op galactische schalen — enorm
- Duizenden melkwegstelsels gemeten, met consistent patroon

---

## 3. Hubble-spanning

### 3.1 De waarneming

De Hubble-constante H0 — de huidige expansiesnelheid van het heelal —
wordt op twee fundamenteel verschillende manieren gemeten, met inconsistente
resultaten:

**Vroeg heelal (CMB):**
- Planck (2018): H0 = 67.4 +/- 0.5 km/s/Mpc
- ACT (2025): bevestigt Planck-waarde, ook uit polarisatiedata

**Laat heelal (lokaal):**
- SH0ES (Cepheïden + supernovae): H0 = 73.04 +/- 1.04 km/s/Mpc
- TDCOSMO (gravitationele lenzing, dec 2025): bevestigt lokale waarde
- Freedman et al. (JWST, 2025): H0 = 70.4 +/- 2.1 km/s/Mpc (overlap met beide)

**Spanning:** Het verschil is ~6 km/s/Mpc, ofwel ~9%. De statistische
significantie is >5 sigma wanneer SH0ES vs. Planck vergeleken wordt.

**Status (2026):** De spanning is NIET opgelost. Nieuwe metingen met JWST
en kosmische lenzen bevestigen afwisselend de hoge en de lage waarde.
De spanning wordt door sommigen een "kosmologische crisis" genoemd.

### 3.2 Hoe de ART het verklaart (of niet)

De ART zelf zegt niets over de waarde van H0 — die is een vrije parameter
in de Friedmann-vergelijkingen. De spanning zit in het standaard kosmologisch
model (LambdaCDM) dat de CMB-data extrapoleert naar het huidige tijdperk.

Als de spanning reeel is (geen systematiek), dan klopt een van de aannames van
LambdaCDM niet:
- De donkere-energievergelijking w = -1 (constante Lambda)
- De Friedmann-vergelijkingen met alleen bekende componenten
- De geluidshorison op het moment van recombinatie

Voorgestelde oplossingen: early dark energy, extra relativistische deeltjes,
interacterende donkere materie/energie, modificatie van de gravitatietheorie.
Geen van deze is overtuigend.

### 3.3 Hoe de ORT het mogelijk anders verklaart

De ORT biedt een **conceptueel ander kader** voor de Hubble-spanning:

**ORT-kosmologie (§12.22-12.25):**
- Het heelal als binnenkant van een zwart gat
- H(t) = c^3 / (2G · M_BH(t)) — de Hubble-parameter volgt uit de BH-massa
- Kosmologische roodverschuiving = gravitationele roodverschuiving

**Mogelijke relevantie voor de Hubble-spanning:**

1. **c_local varieert met positie:** In de ORT is c_local lager nabij massa.
   Het lokale heelal (melkwegstelsels, clusters) zit in een andere
   c_local-omgeving dan de CMB-achtergrond. De lokale H0 zou systematisch
   hoger kunnen zijn dan de globale H0 als de lokale c_local lager is dan
   het gemiddelde.

2. **Schiltheorema op kosmologische schaal:** Lokale metingen van H0 worden
   beinvloed door de massa-dichtheid in de directe omgeving. In de ORT
   verlaagt de omringende kosmische structuur (supercluster, filament)
   de c_local isotroop, wat de lokaal gemeten expansie beinvloedt.

3. **Geen Lambda nodig:** Als de versnelde expansie een gevolg is van de
   BH-groeidynamiek (§12.23-12.24), dan is de extrapolatie van de CMB-data
   anders dan in LambdaCDM. Een ander expansieverloop geeft een andere
   voorspelde lokale H0.

**Eerlijkheid:**
- De ORT-kosmologie is **speculatief** (MODEL.md, §12.25 eerlijkheidstabel)
- Er is geen kwantitatieve voorspelling van H0 uit de ORT
- De bewering "c_local verschilt lokaal" vereist een kwantitatieve berekening
  die nog niet bestaat
- De correspondentie r_s = R_H is wiskundig exact, maar de interpretatie
  als moeder-zwart-gat is een hypothese

### 3.4 Meetprecisie

- Planck CMB: 0.7% precisie op H0
- SH0ES: 1.4% precisie
- TDCOSMO (lenzing): ~2-3% precisie
- Het verschil (~9%) is veel groter dan de individuele fouten
- 5+ sigma significantie — dit is geen randfenomeen

---

## 4. Kosmische versnelling / donkere energie

### 4.1 De waarneming

In 1998 toonden twee teams (SCP en High-z SN Search) aan dat verre
type Ia-supernovae lichtzwakker waren dan verwacht in een vertragende
expansie. Conclusie: de expansie van het heelal versnelt.

Dit vereist een component met negatieve druk: donkere energie, ofwel een
kosmologische constante Lambda met dichtheid:

    rho_Lambda ~ 7 · 10^-30 g/cm^3
    Omega_Lambda ~ 0.685

DESI (Dark Energy Spectroscopic Instrument) heeft in 2025 de meest
precieze BAO-metingen ooit gepubliceerd (DR2, maart 2025):

- 14+ miljoen sterrenstelsels en quasars
- Statistische precisie van 0.24% op de BAO-schaal
- Bij z > 2: 0.65% precisie — beste meting ooit

**Cruciale bevinding:** De DESI DR2-data suggereren dat donkere energie
NIET constant is. De vergelijking w(z) = w0 + wa · (1-a) geeft:
- w0 > -1 (iets minder negatief dan Lambda)
- wa < 0 (w was negatiever in het verleden)

De voorkeur voor dynamische donkere energie is in DR2 niet afgenomen
ten opzichte van DR1 — eerder toegenomen. Maar het is nog niet op
"ontdekkingsniveau" (5 sigma), en er zijn spanningen tussen de BAO-,
CMB- en supernovadatasets onderling.

### 4.2 Hoe de ART het verklaart

De ART staat een kosmologische constante Lambda toe als parameter in de
Einstein-vergelijkingen. Dit is wiskundig consistent maar fysisch
onverklaard:

- **Het finetuning-probleem:** De waargenomen Lambda is ~10^120 keer
  kleiner dan de kwantumveldtheoretische verwachting
- **Het coincidentie-probleem:** Omega_Lambda ~ Omega_m juist in ons
  tijdperk — een toevallige gelijkheid?
- **Dynamische donkere energie** (DESI hint): als w varieert, is Lambda
  niet fundamenteel en is een ander mechanisme nodig

### 4.3 Hoe de ORT het mogelijk anders verklaart

De ORT biedt het meest uitgewerkte alternatief op dit punt (§12.23-12.25):

**1. Versnelde expansie zonder Lambda (§12.23)**

De drie kosmologische era's corresponderen met het groeigedrag van M_BH:

| Era | M_BH(t) | H(t) | a(t) | c_local-gedrag |
|-----|---------|------|------|----------------|
| Inflatie | ~ const (net gevormd) | ~ const (groot) | exp(Ht) | Daalt: hyperexpansie |
| Straling/materie | proportioneel met t | proportioneel met 1/t | t^(1/2) tot t^(2/3) | Stijgt: herstel |
| Huidige era | ~ const (groei stopt) | ~ const (klein) | exp(Ht) | Daalt: hernieuwde expansie |

Het mechanisme: als M_BH stopt met groeien, wordt H constant, en de expansie
wordt exponentieel. Inflatie en de huidige versnelling zijn **hetzelfde
mechanisme** op verschillende schalen. Er is geen Lambda nodig.

**2. Zelfversterkende expansie (positieve terugkoppeling)**

De expansie is zelfversterkend: c_local = 0 op de eventhorizon, de ruimte
ertussen rekt uit, informatie-uitwisseling vertraagt, en de effectieve horizon
groeit sneller dan informatie kan bijhouden. Dit geeft exponentieel gedrag
zonder externe drijfkracht.

**3. Effectieve Lambda uit BH-massa (§12.25)**

    H(t) = c^3 / (2G · M_BH(t))
    Lambda_eff ~ 3c^6 · Omega_Lambda / (4G^2 · M_BH^2)

Met M_BH ~ 9.2 · 10^52 kg levert dit de juiste orde van grootte voor
Lambda ~ 1.1 · 10^-52 m^-2.

**Aansluiting bij DESI:**
De DESI-hint dat w varieert met roodverschuiving past **kwalitatief** bij de ORT:
- w(z) is niet constant omdat de BH-groeidynamiek niet constant is
- De overgang van vertraging naar versnelling is een natuurlijk gevolg van
  veranderend dM_BH/dt
- De ORT voorspelt geen specifieke w0 en wa, maar het concept van een
  varierende w is natuurlijker in de ORT dan in LambdaCDM

**Eerlijkheid:**
- De ORT-kosmologie is **speculatief** — geen kwantitatieve voorspelling van w(z)
- De correspondentie Lambda_eff ~ juiste orde van grootte is suggestief maar
  niet voldoende
- Het finetuning-probleem wordt vervangen door de vraag: waarom is M_BH
  precies deze waarde? (Wat overigens een minder ernstige vraag is dan het
  120-ordes-van-grootte-probleem)

### 4.4 Meetprecisie

- DESI BAO: 0.24% precisie (DR2)
- Type Ia supernovae: ~7% individueel, <1% statistisch
- CMB (Planck): ~0.5% op kosmologische parameters
- De hint voor dynamische w is ~2-3 sigma — suggestief maar niet beslissend

---

## 5. Zwaartekrachtsgolf-echo's

### 5.1 De waarneming

Na een zwart-gat-fusie voorspelt de ART een ringdown: de gevormde Kerr-black-hole
trilt na en dempt exponentieel (quasi-normale modes, QNM). Sommige theorieën
voorspellen daarnaast **echo's** — herhaalde pulsen na de ringdown — als de
eventhorizon niet een perfect absorberende grens is, maar een reflecterend
oppervlak net buiten r_s.

**Status van waarnemingen (2025-2026):**
- GW250114 (jan 2025): SNR = 80, de scherpste gravitatiegolfmeting ooit
  (vergeleken met SNR = 26 voor GW150914)
- Een 2025 model-onafhankelijke zoekactie (arXiv:2512.24730) vond **geen**
  statistisch significante echo's in GW150914, GW231226, of GW250114
- De Hawking-areaalstelling (oppervlak van de horizon neemt nooit af) is
  bevestigd door recente analyses

**Status (2026):** Er zijn tot op heden **geen** overtuigende detecties van
post-merger echo's. De bovengrens op de reflectiecoëfficiënt van de horizon
wordt steeds strenger.

### 5.2 Hoe de ART het verklaart

In de ART is de eventhorizon een perfect absorberende, eenrichtingsgrens.
Informatie die de horizon passeert kan niet terugkeren. Echo's zouden een
schending van dit principe zijn en wijzen op "nieuwe fysica" nabij de horizon
(fuzzballs, firewalls, quantum gravity effects).

Het niet-detecteren van echo's is **consistent** met de ART.

### 5.3 Hoe de ORT het zou behandelen

In de ORT is de eventhorizon de plek waar c_local = 0. Dit geeft een
vergelijkbaar beeld als de ART:

- Bij r = r_s: v_grav = c, c_local = 0
- Alle snelheid gaat naar de massa: er blijft niets over voor beweging
  door ruimte of tijd
- "Zijn = 0" op de horizon — het concept van bestaan verliest betekenis

**De ORT voorspelt GEEN echo's.** De horizon is in de ORT zelfs
"fundamenteler" dan in de ART: het is niet alleen een eenrichtingsgrens,
maar een punt waar de lokale ruimtetijdsnelheid nul is. Reflectie vereist
c_local > 0, wat per definitie niet geldt op de horizon.

**Verschil met de ART:**
De ORT beschrijft het zwart-gat-interieur via de analytische voortzetting
c_interior = c · sqrt(r_s/r - 1) (§12.23). Dit geeft een andere
interpretatie van wat "achter de horizon" is (hyperexpansie), maar verandert
niets aan de buitenkant: geen echo's.

**Ringdown:**
De ORT reproduceert de QNM-frequenties via gelineariseerde perturbaties
van het c_local-veld (§12.26.8), identiek aan de ART. Dit is bevestigd
door de waarnemingen.

**Eerlijkheid:** De niet-detectie van echo's is consistent met zowel de ART
als de ORT. Dit is geen onderscheidend experiment.

### 5.4 Meetprecisie

- LIGO O4: 10x gevoeliger dan O1
- GW250114: SNR = 80
- Echo-detectie vereist SNR >> 100 voor zwakke echo's
- Volgende generatie detectoren (Einstein Telescope, Cosmic Explorer) nodig
  voor definitieve test

---

## 6. Pioneer-anomalie

### 6.1 De waarneming

De Pioneer 10 en 11 ruimtevaartuigen vertoonden een onverklaarde vertraging
van ~8.74 · 10^-10 m/s^2, gericht naar de Zon. De afwijking was constant over
grote afstanden (20-70 AU) en werd voor het eerst gerapporteerd door Anderson
et al. (1998).

### 6.2 De verklaring

In 2012 is de anomalie **volledig verklaard** door anisotrope thermische
straling van de ruimtevaartuigen. De RTG (radioisotopen-thermische
generator) zond warmte uit naar de achterkant van de schotelantenne;
de schotel reflecteerde deze warmte in de vliegrichting. Daarnaast
straalde de elektronica-unit aan de voorkant van het ruimtevaartuig warmte
uit in de vliegrichting. De resulterende fotonendruk gaf een kleine maar
constante vertraging.

Twee onafhankelijke analyses bevestigden dit:
- Turyshev et al. (2012): gedetailleerd thermisch model
- Rievers & Lammerzahl (2011): onafhankelijke bevestiging

**Status (2026):** De Pioneer-anomalie is **opgelost**. De thermische
verklaring is de wetenschappelijke consensus. Er zijn geen nieuwe
ontwikkelingen die deze conclusie in twijfel trekken.

### 6.3 ORT-perspectief

De ORT hoeft hier niets te verklaren — het effect is niet gravitationeel.
De oplossing is een les in systematische effecten: wat jarenlang als
mogelijk nieuwe fysica werd gezien, bleek uiteindelijk een mundaan
thermisch effect.

**Relevantie voor de ORT:** De Pioneer-anomalie herinnert eraan dat
claims van anomalieen streng getoetst moeten worden aan systematische
fouten voordat ze als bewijs voor nieuwe fysica worden ingezet.

### 6.4 Meetprecisie

- Pioneer Doppler-tracking: ~10^-12 m/s^2 precisie
- De anomalie (8.74 · 10^-10 m/s^2) was ~100x groter dan de meetfout
- Thermische modellering reproduceert de waarneming binnen de onzekerheden

---

## 7. Samenvatting en conclusies

### Overzichtstabel

| Anomalie | ART-status | ORT-perspectief | ORT-specifiek mechanisme |
|----------|-----------|-----------------|-------------------------|
| Flyby-anomalie | Onverklaard | Open vraag | Frame-dragging + v_grav-richting (niet uitgewerkt) |
| Rotatiecurves | Via donkere materie | Kwalitatief (speculatief) | Randeffect + multi-mass c_local |
| Hubble-spanning | Onopgelost (>5 sigma) | Kwalitatief (speculatief) | Lokale c_local-variatie |
| Kosmische versnelling | Via Lambda (finetuning-probleem) | Kwalitatief (speculatief) | BH-groeidynamiek, geen Lambda nodig |
| GW-echo's | Niet gevonden (consistent met ART) | Consistent (geen echo's) | c_local = 0 op horizon |
| Pioneer-anomalie | Opgelost (thermisch) | N.v.t. | Niet gravitationeel |

### Waar de ORT potentieel verschil maakt

1. **Kosmische versnelling** — het sterkste punt. De ORT biedt een conceptueel
   schoner alternatief voor Lambda: de versnelling volgt uit BH-groeidynamiek.
   Het finetuning-probleem (10^120 ordes van grootte) verdwijnt. De DESI-hint
   naar dynamische donkere energie past kwalitatief.

2. **Hubble-spanning** — de ORT biedt een conceptueel kader (lokale vs. globale
   c_local) dat de spanning potentieel opheft. Maar er is geen kwantitatieve
   voorspelling.

3. **Rotatiecurves** — de minst uitgewerkte claim. De randeffect-formule is
   kwalitatief en het is niet aangetoond dat het de juiste profiel-vorm geeft.

4. **Flyby-anomalie** — een interessant testgeval voor de ORT, specifiek voor
   de combinatie van v_grav-richting en frame-dragging. Maar er is geen
   berekening.

### Wat de ORT NIET kan

- Geen kwantitatieve voorspellingen voor kosmologische observables
- Geen specifieke rotatiecurve-fits
- Geen voorspelling van H0
- Geen verklaring voor de structuur van de CMB-fluctuaties
- Geen alternatief voor de volledige LambdaCDM-parameterset

### Conclusie

De ORT biedt een **conceptueel aantrekkelijk alternatief kader** voor
enkele van de grootste open vragen in de kosmologie, met name donkere
energie en de Hubble-spanning. De kracht ligt in de eenvoud: een enkel
mechanisme (c_local-dynamiek) vervangt meerdere ad-hoc parameters
(Lambda, donkere materie). De zwakte is het ontbreken van kwantitatieve
voorspellingen die experimenteel getoetst kunnen worden.

De volgende stap voor de ORT zou zijn:
1. Een kwantitatieve voorspelling van w(z) uit de BH-groeidynamiek
2. Een berekening van de flyby-anomalie in het ORT-kader
3. Een rotatiecurve-fit met de randcorrectie-formule

Zonder deze berekeningen blijft de ORT-kosmologie op het niveau van
kwalitatieve hypothesen — suggestief, maar niet falsifieerbaar.

---

## Bronnen

### Flyby-anomalie
- [Flyby anomaly - Wikipedia](https://en.wikipedia.org/wiki/Flyby_anomaly)
- [Flyby Anomaly in the Variation Principle of General Relativity (2024)](https://arxiv.org/abs/2411.12053)
- [The Puzzle of the Flyby Anomaly - ResearchGate](https://www.researchgate.net/publication/225202012_The_Puzzle_of_the_Flyby_Anomaly)
- [Flyby anomaly and frame-dragging - MNRAS](https://academic.oup.com/mnras/article/489/3/3232/5555577)

### Rotatiecurves
- [Galaxy rotation curve - Wikipedia](https://en.wikipedia.org/wiki/Galaxy_rotation_curve)
- [A new empirical fit to galaxy rotation curves (2025)](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2025.1680387/full)
- [Single Parameter Model for Galaxy Rotation Curves (2026)](https://arxiv.org/pdf/2602.24211)
- [Examining Galaxy Rotation and Dark Matter Alternatives (2025)](https://scisimple.com/en/articles/2025-09-17-examining-galaxy-rotation-and-dark-matter-alternatives--a3d2761)

### Hubble-spanning
- [The Hubble tension - CERN Courier](https://cerncourier.com/a/the-hubble-tension/)
- [New cosmic lens measurements deepen Hubble tension (dec 2025)](https://www.sciencedaily.com/releases/2025/12/251209043036.htm)
- [The Hubble Constant and the Crisis in Cosmology: 2025 Status Report](https://newspaceeconomy.ca/2025/12/06/the-hubble-constant-and-the-crisis-in-cosmology-a-2025-status-report/)
- [Hubble tension still unresolved - Ethan Siegel](https://medium.com/starts-with-a-bang/the-hubble-tension-still-unresolved-despite-new-measurements-413a6621b91d)
- [Hubble Constant and Tension - NASA](https://science.nasa.gov/mission/hubble/science/science-behind-the-discoveries/hubble-constant-and-tension/)

### Kosmische versnelling / donkere energie
- [DESI DR2 Results Guide (maart 2025)](https://www.desi.lbl.gov/2025/03/19/desi-dr2-results-march-19-guide/)
- [Dynamical dark energy in light of DESI DR2 - Nature Astronomy](https://www.nature.com/articles/s41550-025-02669-6)
- [Key Results from DESI DR2 - Astrobites](https://astrobites.org/2025/10/06/desi-dr2-part1/)
- [Did DESI DR2 truly reveal dynamical dark energy? - EPJC](https://link.springer.com/article/10.1140/epjc/s10052-025-15076-y)

### Zwaartekrachtsgolf-echo's
- [Model-independent search of gravitational wave echoes (2025)](https://arxiv.org/abs/2512.24730)
- [GW250114 - Wikipedia](https://en.wikipedia.org/wiki/GW250114)
- [Gravitational wave analysis confirms merging black holes (2025)](https://phys.org/news/2025-09-gravitational-analysis-theory-merging-black.html)
- [Constraining Black Hole Horizon Properties (2025)](https://arxiv.org/html/2511.06536)

### Pioneer-anomalie
- [Pioneer Anomaly Solved! - The Planetary Society](https://www.planetary.org/articles/3459)
- [Support for the Thermal Origin - PRL (2012)](https://link.aps.org/doi/10.1103/PhysRevLett.108.241101)
- [Pioneer anomaly - Wikipedia](https://en.wikipedia.org/wiki/Pioneer_anomaly)
