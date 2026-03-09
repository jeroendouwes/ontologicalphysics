# Overwegingen: Isotropie van ruimtelijke uitrekking in de ORT

*Werkdocument — niet voor publicatie*

---

## 1. De kernvraag

De ORT stelt dat c_local een **scalair** is: het snelheidsbudget daalt nabij massa in alle richtingen gelijk. Tijddilatatie is inderdaad isotroop — een klok tikt langzamer ongeacht zijn oriëntatie. Maar in de ART (Schwarzschild) is de ruimtelijke uitrekking **anisotoop**: alleen radiaal.

**Waarom zou ruimte zich anders gedragen dan tijd?**

---

## 2. Wat de wiskunde zegt

### ART (Schwarzschild-coördinaten)
- g_tt = (1 - r_s/r) — isotroop (scalair)
- g_rr = 1/(1 - r_s/r) — radiaal uitgerekt
- g_θθ = r² — tangentieel NIET uitgerekt
- Omtrek = 2πr (onveranderd)

### ORT als isotroop geïnterpreteerd
- g_tt = (c_local/c)² = (1 - r_s/r) — zelfde
- g_rr = (c/c_local)² = 1/(1 - r_s/r) — zelfde
- g_θθ = r²/(1 - r_s/r) — WEL uitgerekt
- Omtrek = 2πr/√(1 - r_s/r) → divergeert bij r_s

### Consequenties van isotrope uitrekking
- Fotonsfeer verschuift: 1.5 r_s (ART) → 2.0 r_s (ORT isotroop)
- Zwart-gatschaduw: b_crit = 2.60 r_s (ART) → 4.0 r_s (ORT) — **54% groter**
- Tangentiële lichtsnelheid: c√(1-r_s/r) (ART) → c(1-r_s/r) (ORT) — **langzamer**

---

## 3. EHT-waarnemingen: bewijs of bevestigingsbias?

### Het standaardverhaal
De EHT meet de schaduw van M87* als ~42 μas, "consistent met ART". De onzekerheid is ~10-17%. Een 54% afwijking zou buiten deze marge vallen → isotrope ORT uitgesloten?

### Kanttekeningen

**a) Modelafhankelijke reconstructie**
De EHT maakt geen directe foto. Het reconstrueert een beeld uit sparse interferometriedata (slechts ~20 telescopen, onvolledig u-v vlak). De reconstructie gebruikt:
- GRMHD-simulaties gebaseerd op de Kerr-metriek (ART)
- Ray-tracing met ART-geodeten
- Bayesiaanse modelfit met ART als prior

Als je de data zou reconstrueren met een ORT-model (andere geodeten, andere lensing), zou het gereconstrueerde beeld anders zijn.

**b) Massa-schatting is metriekafhankelijk**
De massa van M87* (6.5 × 10⁹ M☉) komt uit sterrenbewegingen en gasdynamica, geanalyseerd met Newtoniaanse/ART-modellen. In de ORT zou dezelfde stellaire kinematica een andere massa-schatting opleveren. Een lagere massa → kleinere r_s → kleinere verwachte schaduw, wat het 54%-verschil deels zou compenseren.

**c) Wetenschappelijke bevestigingsbias**
Wetenschappers neigen ernaar te zien wat ze verwachten. De hele EHT-pipeline is ontworpen en gevalideerd tegen ART-voorspellingen. Een systematische bias is niet uitgesloten.

**d) Resolutiebeperking**
De EHT heeft een resolutie van ~20 μas. De schaduw van ~42 μas is slechts ~2 resolutie-elementen breed. De precieze schaduwgrootte is minder zeker dan de publicaties suggereren.

---

## 4. Het argument voor anisotrope uitrekking

### De bron-richting is speciaal
Zwaartekracht is een **radiaal** veld. De ontsnappingssnelheid v_grav wijst naar de massa. Er IS een voorkeursrichting in de ruimte: de lijn bron→waarnemer.

- Tijd heeft geen voorkeursrichting → tijddilatatie is isotroop
- Ruimte heeft een voorkeursrichting (radiaal) → uitrekking is anisotoop

### Lading bevestigt het patroon
In Reissner-Nordström (geladen massa) verandert f(r) = 1 - r_s/r + r_Q²/r². Maar ook hier: alleen g_tt en g_rr worden beïnvloed, g_θθ = r² blijft onveranderd. Lading is ook radiaal → dezelfde anisotropie.

### v_grav als vector
c_local is scalair (het resterende snelheidsbudget), maar v_grav is een **vector** (gericht naar de massa). De ruimtelijke uitrekking wordt veroorzaakt door v_grav, niet door c_local zelf. En een vector heeft een richting.

Interpretatie: c_local bepaalt de tijddilatatie (isotroop, scalair). Maar de ruimtelijke uitrekking komt uit de **gradiënt** van c_local, die radiaal is.

---

## 5. Het argument voor isotrope uitrekking

### c_local is fundamenteel
In de ORT is c_local het lokale snelheidsbudget. Een liniaal in elke richting wordt "gemeten" met lokaal licht. Als licht in alle richtingen met c_local beweegt, dan is elke lokale meting — inclusief afstanden — in alle richtingen gelijk beïnvloed.

### Isotrope coördinaten bevestigen het
De Schwarzschild-metriek in isotrope coördinaten (Eddington) IS conform vlak:
dl² = (1 + r_s/(4ρ))⁴ (dρ² + ρ²dΩ²)

Dit is dezelfde factor in alle richtingen. De "anisotropie" in Schwarzschild-coördinaten is een coördinaat-artefact.

### Binnen het zwart gat: isotropie hersteld
Binnen de eventhorizon (r < r_s) wisselen r en t van karakter. Er is geen radiale "bron-richting" meer. De expansie is isotroop — precies het kosmologisch principe (CMB isotropie, homogene expansie).

Als het heelal de binnenkant van een zwart gat is (ORT-hypothese), dan is isotropie de NATUURLIJKE toestand. De anisotropie buiten het zwart gat is het speciale geval, niet andersom.

### Licht draait niet om het heelal
Als de ruimte binnen het zwart gat isotroop uitrekt, kan licht niet om het heelal heen draaien — de expansie "wint" van de lichtsnelheid. Dit is consistent met de waarneming dat we geen herhaalde beelden van het heelal zien.

---

## 6. Mogelijke synthese

Misschien is het geen of/of maar een kwestie van perspectief:

| | Buiten BH | Binnen BH |
|---|---|---|
| Bron-geometrie | Radiaal (massa op r=0) | Geen centrum |
| v_grav | Vector (radiaal) | Geen voorkeursrichting |
| Ruimtelijke uitrekking | Anisotoop (radiaal) | Isotroop |
| Metriek | Schwarzschild | FLRW-achtig |

c_local is scalair. Maar de **manier waarop** c_local de ruimte beïnvloedt hangt af van de brongeometrie:
- Puntbron → radiaal veld → anisotrope uitrekking
- Geen bron (binnenkant) → geen voorkeursrichting → isotrope uitrekking

Dit zou consistent zijn met:
- ART buiten BH (Schwarzschild: anisotoop) ✓
- ART binnen BH (FLRW: isotroop) ✓
- EHT-waarnemingen (fotonsfeer bij 1.5 r_s) ✓
- Kosmologisch principe (isotrope expansie) ✓

---

## 7. Wat moet er in §4.6 staan?

### Optie A: Isotroop claim handhaven
- Eerlijk presenteren: ORT voorspelt isotrope uitrekking
- Consequenties benoemen: andere fotonsfeer, grotere schaduw
- EHT-data bespreken als potentieel tegenbewijs, met kanttekeningen over modelafhankelijkheid
- Sterkste versie van de theorie: laat de data beslissen

### Optie B: Anisotoop accepteren, isotropie laten vallen
- ORT reproduceert Schwarzschild exact (g_tt en g_rr)
- g_θθ = r² accepteren (niet isotroop uitgerekt)
- Verklaring: v_grav is een vector → ruimtelijke anisotropie
- Binnen BH: isotropie hersteld (kosmologie-notebook)

### Optie C: Nuance — coördinaatafhankelijk
- In isotrope coördinaten IS de Schwarzschild-ruimte isotroop
- In Schwarzschild-coördinaten niet
- De ORT-claim is correct in isotrope coördinaten
- Geen fysisch verschil met ART — het IS dezelfde metriek
- Maar: dan is er ook geen "testbare voorspelling" bij 2e orde

---

## 8. Open vragen voor morgen

1. Is de 2e-orde verschil (3/8 vs 5/16) een echt fysisch verschil of een coördinaat-artefact?
2. Kan de EHT-data opnieuw geanalyseerd worden met een ORT-metriek?
3. Hoe verhoudt de fotonsfeer bij 2r_s zich tot Sgr A* (waar de massa nauwkeuriger bekend is)?
4. Als v_grav de anisotropie verklaart: hoe formaliseer je dat wiskundig?
5. Hoe gaat de buitenkant (anisotoop) over in de binnenkant (isotroop) bij r = r_s?
6. Licht kan niet om het heelal draaien — is dit een consequentie van isotrope uitrekking of van de expansiesnelheid?
