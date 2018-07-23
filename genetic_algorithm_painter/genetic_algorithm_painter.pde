// Processing 3.3.7.
// 23.07.2018
// GitHub repo creation and refactoring and translation to English
// 19.07.2018 avaruustursas
// 21.7.2018
//
// 21.7.2018
// One gene is following tuple (x1,y1,x2,y2,x3,y3,r,g,b,a)
// Values for xi,yi,r,g,b,a are in range [0.0,1.0] in the gene
// Fitness is measured pixel by pixel comparison
// Euclidian distance between colours; (255,255,255) is normalized to (1,1,1)
// Distance of every pixel from the input to current candidate is calculated.
// If the result is lesser than current best (or shortest distance) then the
// best genome is updated to the new one
//
// Generation: mutate the fittest result 50 times and update the best if necessary


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


FloatList[] mutate_genome(FloatList[] genome) {
  // first test was to just .shuffle the genome.
  // The results were too random
  FloatList[] uusi_genomi = new FloatList[NGENES];
  uusi_genomi = copy_genome(genome);
  for (FloatList g : uusi_genomi) {
    float rand = random(0.0,1.0);
    if (rand > 0.99)
      g.shuffle();
    float ehto = random(0.0,1.0);
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
  //print("Starting to compare image vs screen\n");
  img.loadPixels();
  color[] img_pixels = img.pixels;
  color[] screen_pixels = load_generated_pixels_from_screen(256+256+2,0,256,256); // screen
  float sum_value = 0;
  for (int i=0;i<img_pixels.length;i++) {
    float delta_r = (img_pixels[i] >> 16 & 0xFF)-(screen_pixels[i] >> 16 & 0xFF);
    float delta_g = (img_pixels[i] >> 8 & 0xFF)-(screen_pixels[i] >> 8 & 0xFF  );
    float delta_b = (img_pixels[i] & 0xFF)-(screen_pixels[i] & 0xFF);
    sum_value += sqrt((sq(delta_r)+sq(delta_g)+sq(delta_b))/(255*255))/sqrt(3); // sqrt and refactoring for norm 23.7.2018 11:18
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
int generation = 0;
int generation_of_the_fittest = 0;

void setup() {
  size(800,600);
  background(0);
  // load input image (256x256 pixels)
  img = loadImage("testikuva.png");
  //img = loadImage("lena.png");
  genome = create_genome(NGENES);
  print("setup ready\n");
}

FloatList[] fittest_one = new FloatList[NGENES];
float fittest_value = 1.0;
IntList steps_to_new_one = new IntList();

void draw() {
  background(0);
  image(img,0,0);
  // Generate 50 new mutated genomes from the fittest one
  for(int u=0;u<50;u++) {
    background(0);
    image(img,0,0);
    FloatList[] mutaatio = mutate_genome(genome);
    if(random(0,100) > 96) {
      mutaatio = switch_gene_places(mutaatio);
    }
    // draw the mutated version to get something to look at
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
    if (evaluation < fittest_value) {
      fittest_value = evaluation;
      fittest_one = copy_genome(mutaatio);
      steps_to_new_one.append(generation-generation_of_the_fittest);
      generation_of_the_fittest = generation;
    }
  }
  // Draw the fittest version
  for (FloatList gene : fittest_one) {
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
  print("------\n");
  // update the fittest genome for seed of the next generation
  genome = copy_genome(fittest_one);
  print("Ongoing generation: "+generation+", generation of the fittest ("+generation_of_the_fittest+")\n");
  print("Average amount of steps to generate better fit "+karvo(steps_to_new_one)+" (median "+mediaani(steps_to_new_one)+") generations\n");
  generation++;
}
