#include "LU.h"
float A_random[N][N]={-0.0363885238766670227050781250000000000000000000000000,
-0.0320882722735404968261718750000000000000000000000000,
-0.0014431107556447386741638183593750000000000000000000,
-0.0876644775271415710449218750000000000000000000000000,
-0.0509621202945709228515625000000000000000000000000000,
0.0322982817888259887695312500000000000000000000000000,
0.0578212440013885498046875000000000000000000000000000,
-0.1164292618632316589355468750000000000000000000000000,
-0.1506015211343765258789062500000000000000000000000000,
0.2177740335464477539062500000000000000000000000000000,
-0.2168580442667007446289062500000000000000000000000000,
-0.1760680377483367919921875000000000000000000000000000,
0.1623535752296447753906250000000000000000000000000000,
-0.1465168148279190063476562500000000000000000000000000,
0.4816236197948455810546875000000000000000000000000000,
-0.0412726998329162597656250000000000000000000000000000,
0.0874162167310714721679687500000000000000000000000000,
0.0617984458804130554199218750000000000000000000000000,
-0.0011459626257419586181640625000000000000000000000000,
-0.0886914879083633422851562500000000000000000000000000,
0.0589736737310886383056640625000000000000000000000000,
-0.0042359791696071624755859375000000000000000000000000,
0.1591832935810089111328125000000000000000000000000000,
-0.0505535565316677093505859375000000000000000000000000,
-0.0633621737360954284667968750000000000000000000000000,
-0.0487045422196388244628906250000000000000000000000000,
0.0488059930503368377685546875000000000000000000000000,
0.0573621056973934173583984375000000000000000000000000,
-0.0162379760295152664184570312500000000000000000000000,
0.0450068302452564239501953125000000000000000000000000,
-0.0865892171859741210937500000000000000000000000000000,
0.0216104816645383834838867187500000000000000000000000,
-0.0732838064432144165039062500000000000000000000000000,
0.1196968257427215576171875000000000000000000000000000,
-0.1047580689191818237304687500000000000000000000000000,
-0.0676544457674026489257812500000000000000000000000000,
0.0183143056929111480712890625000000000000000000000000,
-0.0416854396462440490722656250000000000000000000000000,
0.2704969644546508789062500000000000000000000000000000,
-0.0738006979227066040039062500000000000000000000000000,
0.1950599998235702514648437500000000000000000000000000,
-0.0242610014975070953369140625000000000000000000000000,
0.0964509546756744384765625000000000000000000000000000,
0.0406045056879520416259765625000000000000000000000000,
0.0915637612342834472656250000000000000000000000000000,
0.0232575405389070510864257812500000000000000000000000,
-0.1068224906921386718750000000000000000000000000000000,
0.0084624132141470909118652343750000000000000000000000,
-0.0142460400238633155822753906250000000000000000000000,
0.1363401710987091064453125000000000000000000000000000,
-0.1216815337538719177246093750000000000000000000000000,
0.1230439767241477966308593750000000000000000000000000,
0.0920307785272598266601562500000000000000000000000000,
-0.1211169883608818054199218750000000000000000000000000,
-0.0215412769466638565063476562500000000000000000000000,
0.1837451159954071044921875000000000000000000000000000,
0.1248173639178276062011718750000000000000000000000000,
0.0406893864274024963378906250000000000000000000000000,
0.2133611887693405151367187500000000000000000000000000,
0.4022989571094512939453125000000000000000000000000000,
0.1926660835742950439453125000000000000000000000000000,
0.0246231257915496826171875000000000000000000000000000,
-0.4380384087562561035156250000000000000000000000000000,
0.2951385378837585449218750000000000000000000000000000,
};
float b_random[N]={-0.0958733931183815002441406250000000000000000000000000,
-0.0291961226612329483032226562500000000000000000000000,
-0.0776031389832496643066406250000000000000000000000000,
0.0723175778985023498535156250000000000000000000000000,
0.1002933457493782043457031250000000000000000000000000,
0.2361856848001480102539062500000000000000000000000000,
-0.2703775465488433837890625000000000000000000000000000,
0.2346844077110290527343750000000000000000000000000000,
};
