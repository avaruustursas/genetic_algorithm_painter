// 19.07.2018 avaruustursas
// 21.7.2018
// Tehdään 800,600 ikkuna
// 256*256 kuva ladataan vasempaan yläkulmaan
// sen alle paras tulos
// vieressä vaikka 2x2 eri versiota
//
// 21.7.2018
// geeni (x1,y1,x2,y2,x3,y3,r,g,b,a)
// nuo arvoit voivat olla myös aina 0.0-1.0
// vertailu
// matriisit eri värikanavien erotuksista
// normalisointi?
// sitten kerrotaan tulos vektorillla (1,1,1,1,..,1) -> vektori
// sitten vektorin komponentit yhteen?
//
// TODO:
// genomitaulukko
// kuvan vertailu
// mutaatio

int NGENES = 150;

// generate one gene
FloatList create_gene() {
  int n = 10;
  FloatList new_gene = new FloatList();
  for(int i=0;i<n;i++)
    new_gene.append(random(0.0,1.0));
  return new_gene;
}

FloatList[] create_genome(int n) {
  print("Init new random genome (len="+n+")\n");
  FloatList[] genome = new FloatList[n];
  for(int i=0;i<genome.length;i++)
    genome[i] = create_gene();
  return genome;
}

FloatList[] copy_genome(FloatList[] vanha) {
  FloatList[] uusi = new FloatList[NGENES];
  int ind = 0;
  for (FloatList g : vanha) {
    FloatList uusi_g = new FloatList();
    for (float ug : g) {
      uusi_g.append(ug);
    }
    uusi[ind] = uusi_g;
    ind++;
  }
  return uusi;
}

//TODO
// lisää jokin tn, että pari kolmiota vaihtaa järjestystä genome listassa

FloatList[] mutate_genome(FloatList[] genome) {
  // eka naiivi, FloatList.shuffle
  // oli liian random
  FloatList[] uusi_genomi = new FloatList[NGENES];
  uusi_genomi = copy_genome(genome);
  for (FloatList g : uusi_genomi) {
    float rand = random(0.0,1.0);
    if (rand > 0.99)
      g.shuffle();
    float ehto = random(0.0,1.0);
    //if (random(0.0,1.0) < 0.2) {
    if (ehto < 0.2) {
      // random valid g index
      int ind = constrain(int(random(0,g.size())),0,g.size()-1);
      float old_value = g.get(ind);
      float sign = random(-1,1);
      //g.set(ind,constrain(old_value+sign*random(0.0,0.5),0.0,1.0));
      if (sign <= 0) {
        g.set(ind,constrain(old_value*0.9,0.0,1.0));
      } else {
        g.set(ind,constrain(old_value*1.1,0.0,1.0));
      }
    }
    if (random(0.0,1.0) <= 0.02) {
      // random valid g index
      int ind = constrain(int(random(0,g.size())),0,g.size()-1);
      float old_value = g.get(ind);
      float sign = random(-1,1);
      if (sign <= 0) {
        float new_value = old_value-old_value/2.0;
        g.set(ind,constrain(new_value,0.0,1.0));
        //g.set(ind,constrain(old_value+sign*random(0.0,0.5),0.0,1.0)/2.0);
      }
      else {
        float new_value = old_value+old_value/2.0;
        g.set(ind,constrain(new_value,0.0,1.0));
      }
    }
    if (random(0.0,1.0) <= 0.01) {
      // random valid g index
      int ind = constrain(int(random(0,g.size())),0,g.size()-1);
      float old_value = g.get(ind);
      float sign = random(-1,1);
      if (sign <= 0) {
        float new_value = old_value-1.0;
        g.set(ind,constrain(new_value,0.0,1.0));
        //g.set(ind,constrain(old_value+sign*random(0.0,0.5),0.0,1.0)/2.0);
      }
      else {
        float new_value = old_value+1.0;
        g.set(ind,constrain(new_value,0.0,1.0));
      }
    }
  }
  return uusi_genomi;
}

FloatList[] switch_gene_places(FloatList[] genome) {
  int l = genome.length;
  int ind1 = int(random(0,l-1));
  int ind2 = int(random(0,l-1));
  FloatList temp = genome[ind1];
  genome[ind1] = genome[ind2];
  genome[ind2] = temp;
  return genome;
}

color[] load_generated_pixels_from_screen(int x0, int y0, int w, int h) {
  //print("Loading pixels from screen\n");
  loadPixels();
  int ind = 0;
  color[] screen_pixels = new color[w*h];
  for (int i=y0;i<y0+h;i++) {
    for(int j=x0;j<x0+w;j++) {
      screen_pixels[ind] = get(j,i);
      ind++;
    }
  }
  //print("Done with loading pixels from screen\n");
  return screen_pixels;
}

float compare_genome_to_image(PImage img) {
  //float result;
  //print("Starting to compare image vs screen\n");
  img.loadPixels();
  color[] img_pixels = img.pixels;
  color[] screen_pixels = load_generated_pixels_from_screen(256+256+2,0,256,256); // screen
  float sum_value = 0; // todella naiivi
  for (int i=0;i<img_pixels.length;i++) {
    //float delta_r = red(img_pixels[i])-red(screen_pixels[i]);
    //float delta_g = green(img_pixels[i])-green(screen_pixels[i]);
    //float delta_b = blue(img_pixels[i])-blue(screen_pixels[i]);
    float delta_r = (img_pixels[i] >> 16 & 0xFF)-(screen_pixels[i] >> 16 & 0xFF);
    float delta_g = (img_pixels[i] >> 8 & 0xFF)-(screen_pixels[i] >> 8 & 0xFF  );
    float delta_b = (img_pixels[i] & 0xFF)-(screen_pixels[i] & 0xFF);
    // lisää vielä neliöjuuri niin saadaan oikeasti vektorin pituus 22.7.2018 02:26
    sum_value += sqrt((sq(delta_r)+sq(delta_g)+sq(delta_b))/(255*255))/sqrt(3); // sqrt and refactoring for norm 23.7.2018 11:18
    //sum_value += (abs(delta_r)+abs(delta_g)+abs(delta_b))/(3*255); // 3 added 22.7.2018 02:13
  }
  //print("Comparison done\n");
  return sum_value/(256*256);
}

float karvo(IntList l) {
  int summa = 0;
  for(int i : l) {
    summa += i;
  }
  return float(summa)/l.size();
}

float mediaani(IntList l) {
  l.sort();
  if (l.size() % 2 == 0) {
    int ind = l.size()/2;
    int ind1 = ind-1;
    int ind2 = ind+1;
    float value = (l.get(ind1)+l.get(ind2))/2.0;
    return value;
  }
  else {
    return l.get(floor(l.size()/2));
  }
}

PImage img;
FloatList[] genome;
int sukupolvi = 0;
int parhaan_sukupolvi = 0;
void setup() {
  size(800,600);
  background(0);
  // load input image (256x256 pixels)
  //img = loadImage("testikuva.png");
  img = loadImage("lena.png");
  //print(genome);
  //print("\n");
  //print("Init new random genome\n");
  // ei toimi
  //for(FloatList g : genome)
  //  g = create_gene();
  //for(int i=0;i<genome.length;i++)
  //  genome[i] = create_gene();
  //print(genome);
  //print("\n");
  genome = create_genome(NGENES);
  print("setup valmis\n");
}
FloatList[] paras = new FloatList[NGENES];
float paras_value = 1.0;
IntList parempi_askeleet = new IntList();
void draw() {
  // lataa kuva ja pikselit matriisiin
  // ja laske erotus matriisioperaatioilla
  //fill(123,123,223,1.0);
  //triangle(400,100,700,500,100,500);
  //float r = random(0,255);
  //float g = random(0,255);
  //float b = random(0,255);
  //float alpha = random(0,100);
  //float x1 = random(0,800);
  //float x2 = random(0,800);
  //float x3 = random(0,800);
  //float y1 = random(0,600);
  //float y2 = random(0,600);
  //float y3 = random(0,600);
  
  background(0);
  image(img,0,0);
  // tee 50 eri genomia ja valitse paras niistä, mikä piirretään
  //FloatList[] paras = new FloatList[100];
  //float paras_value = 1.0;
  for(int u=0;u<50;u++) {
    background(0);
    image(img,0,0);
    //FloatList[] genome = create_genome(100);
    FloatList[] mutaatio = mutate_genome(genome);
    if(random(0,100) > 96) {
      mutaatio = switch_gene_places(mutaatio);
    }
    for (FloatList gene : mutaatio) {
    noStroke();
    float x1 = (width-(256+256+2))*gene.get(0)+256+256+2;
    float y1 = (256)*gene.get(1);
    float x2 = (width-(256+256+2))*gene.get(2)+256+256+2;
    float y2 = (256)*gene.get(3);
    float x3 = (width-(256+256+2))*gene.get(4)+256+256+2;
    float y3 = (256)*gene.get(5);
  
    float r = 255*gene.get(0);
    float g = 255*gene.get(1);
    float b = 255*gene.get(2);
    float alpha = 255*gene.get(3);
    fill(r,g,b,alpha);
    
    triangle(x1,y1,x2,y2,x3,y3);
    }
    float evaluation = compare_genome_to_image(img);
    if (evaluation < paras_value) {
      paras_value = evaluation;
      paras = copy_genome(mutaatio);
      // entinen parhaan ja uuden sukupolvien erotus talteen
      parempi_askeleet.append(sukupolvi-parhaan_sukupolvi);
      parhaan_sukupolvi = sukupolvi;
      //genome = paras;
      //print("uusi paras löytynyt\n");
    }
  }
  //FloatList[] genome = create_genome(100);
  //float evaluation = compare_genome_to_image(img);
  // second quarter
  //print("Piirretään paras\n");
  for (FloatList gene : paras) {
    noStroke();
    float x1 = (width-(256+256+2))*gene.get(0)+256+256+2;
    float y1 = (height-256)+gene.get(1)*256;
    float x2 = (width-(256+256+2))*gene.get(2)+256+256+2;
    float y2 = (height-256)+gene.get(3)*256;
    float x3 = (width-(256+256+2))*gene.get(4)+256+256+2;
    float y3 = (height-256)+gene.get(5)*256;
  
    float r = 255*gene.get(0);
    float g = 255*gene.get(1);
    float b = 255*gene.get(2);
    float alpha = 255*gene.get(3);
    fill(r,g,b,alpha);
    
    triangle(x1,y1,x2,y2,x3,y3);
  }
  //print("genomi valmis, aloitetaan evaluaatio\n");
  //float evaluation = compare_genome_to_image(img);
  //print(evaluation);
  //print("\n");
  // full screen
  //for (FloatList gene : genome) {
  //  noStroke();
  //  float x1 = width*gene.get(0);
  //  float y1 = height*gene.get(1);
  //  float x2 = width*gene.get(2);
  //  float y2 = height*gene.get(3);
  //  float x3 = width*gene.get(4);
  //  float y3 = height*gene.get(5);
  
  //  float r = 255*gene.get(0);
  //  float g = 255*gene.get(1);
  //  float b = 255*gene.get(2);
  //  float alpha = 255*gene.get(3);
  //  fill(r,g,b,alpha);
    
  //  triangle(x1,y1,x2,y2,x3,y3);
  //}
  
  //delay(10);
  print("------\n");
  genome = copy_genome(paras);
  print("Sukupolvi: "+sukupolvi+", paras ("+parhaan_sukupolvi+")\n");
  print("Parempi löytyy keskimäärin "+karvo(parempi_askeleet)+" (mediaani "+mediaani(parempi_askeleet)+") sukupolvessa\n");
  sukupolvi++;
}
