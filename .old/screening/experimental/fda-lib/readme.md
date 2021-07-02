# FDA lib screening 

## aim
Screen FDA lib against BM3 mutants - as many as possible!

## planning

I transferred the SelleckChem FDA lib from vials to PCR plates to avoid de-lidding again (preserving layout). I have roughly 150-200 µl of each (~950).

The first mutants I want to screen are:
1) A82F/F87V
2) A82F 
3) wt

Then whichever mutants I can get made after that. Sam and David say that I should save my money on the MolPort library and instead spend it on getting sequences synthesized, which will make this whole process much faster.

From previous runs with 96 compounds, it takes about 30 mins to dispense compounds into a plate with the echo and maybe 15-20 to dispense the DMSO. 950/48 = 20 plates per mutant, plus the blanks. That's 10h for compounds and 5h for DMSO. If I use the Cai group echo then I'll need to be there to replace the plates each time they finish, so splitting it over 2 days is essential.

To make this viable in the long run, I'll need to use the FBRH platform for dispensing. I think Mark Dunstan is the guy to talk to. I'll do the double mutant manually to proove the concept.

I need to make some changes to ```echo``` to make sure that it works ok.

The screening process itself can be done with manual multichannel pipetting, though the calibration of each channel will become an issue. Using injection in the BMG platereader is also an option, the errors will be constant at least. In the past, bubble formation with the injector hasn't been an issue, but we'll see - now the buffer contains BSA as a stablizer which can form bubbles sometimes, in which case, the plates will need to be spun. In an ideal world, I'd do the entire screen on the FBRH platform, which I think has a build in plate centrifuge.
