#!/usr/bin/ruby -w

# Open the file with the AMBER pots
file = File.new("all_amber_pots.txt", "r")

# Open the file with the GULP equivalence
# This file consists of lines containing AMBER FF type = GULP type
equiv_file = File.new("equiv.txt", "r")

equiv_table = Hash.new("");
while (line=equiv_file.gets)
  tokens = line.split
  equiv_table[tokens[0]] = tokens[2]
end

equiv_table.each do|fftype,symbol|
  puts "# %2s =  %2s" % [fftype, symbol]
end
puts

atom_labels = Hash.new(0)
atom_epsilon = Hash.new("");
atom_sigma = Hash.new("");

bond_labels = Hash.new(0);
bond_k = Hash.new("");
bond_r = Hash.new("");

three_labels = Hash.new(0);
three_k = Hash.new("");
three_angle = Hash.new("");

tors_labels = Hash.new(0);
tors_k = Hash.new("");
tors_phase = Hash.new("");
tors_n = Hash.new("");

improper_labels = Hash.new(0);
improper_k = Hash.new("");
improper_phase = Hash.new("");
improper_n = Hash.new("");

while (line = file.gets)
  line.strip!
  if line == "ATOM    RES  RESNAME  NAME  TYPE   LJ Radius    LJ Depth      Mass    Charge GB Radius GB Screen"
    have_atom = true
  elsif line == "Atom 1               Atom 2               R eq       Frc Cnst"
    have_bond = true
  elsif line == "Atom 1               Atom 2               Atom 3               Frc Cnst   Theta eq"
    have_angle = true
  elsif line == "Atom 1               Atom 2               Atom 3               Atom 4                Height     Periodic.  Phase      EEL Scale  VDW Scale"
    have_tors = true
  else 
    # change brackets to spaces and then can split on whitespace
    line.gsub! "(", " "
    line.gsub! ")", " "
#   puts line
    tokens = line.split
    if (have_atom)
      if tokens[4].nil?
        have_atom = false;
      else
        atom = tokens[4]
        atom_labels[atom]+=1
        atom_sigma[atom] = tokens[5]
        atom_epsilon[atom] = tokens[6]
      end
    elsif (have_bond)
      if tokens[2].nil?
        have_bond = false;
      else
        if (tokens[2] <=> tokens[5]) < 0
          pair = tokens[2] + " " + tokens[5]
        else
          pair = tokens[5] + " " + tokens[2]
        end
        bond_labels[pair]+=1
        bond_k[pair] = tokens[7];
        bond_r[pair] = tokens[6];
      end
    elsif (have_angle)
      if tokens[2].nil?
        have_angle = false
      else
        if (tokens[2] <=> tokens[8]) < 0
          triple = tokens[2] + " " + tokens[5] + " " + tokens[8]
        else
          triple = tokens[8] + " " + tokens[5] + " " + tokens[2]
        end
        three_labels[triple]+=1
        three_k[triple] = tokens[9];
        three_angle[triple] = tokens[10];
      end
    elsif (have_tors)
      if tokens[0].nil?
        have_tors = false
      else
        if tokens[0] =~ /[[:alpha:]]/
          if tokens[0] == "I"
            have_improper = true
          end
          # TODO - do I have to do anything with other codes?
          # so remove first character of line and retokenise
          line.slice!(0)
          tokens = line.split
        end
        if (tokens[2] <=> tokens[11]) < 0
          tuple = tokens[2] + " " + tokens[5] + " " + tokens[8] + " " + tokens[11]
        else
          tuple = tokens[11] + " " + tokens[8] + " " + tokens[5] + " " + tokens[2]
        end
        if ((tokens[2] <=> tokens[11]) == 0) && ((tokens[5] <=> tokens[8]) > 0)
          tuple = tokens[2] + " " + tokens[5] + " " + tokens[8] + " " + tokens[11]
        end
        tuple = tuple + " " + tokens[13]
        if have_improper
          improper_labels[tuple]+=1
          if (tokens[14].to_f > 180.2) && (tokens[14].to_f < 179.8) 
            puts "WARNING: improper torsion angle is not 180 degrees"
          end

          if improper_k[tuple] == ""
            improper_k[tuple] = tokens[12]
            improper_phase[tuple] = tokens[14]
            improper_n[tuple] = tokens[13]
          elsif (improper_k[tuple] == tokens[12]) && (improper_phase[tuple] == tokens[14]) && (improper_n[tuple] == tokens[13])
          else 
            puts "WARNING: %9s %8s %8s %8s" %  [tuple, improper_k[tuple], improper_phase[tuple], improper_n[tuple]]
            improper_k[tuple] = tokens[12]
            improper_phase[tuple] = tokens[14]
            improper_n[tuple] = tokens[13]
          end
          have_improper = false
        else
          tors_labels[tuple]+=1
          if tors_k[tuple] == ""
            tors_k[tuple] = tokens[12]
            tors_phase[tuple] = tokens[14]
            tors_n[tuple] = tokens[13]
          elsif (tors_k[tuple] == tokens[12]) && (tors_phase[tuple] == tokens[14]) && (tors_n[tuple] == tokens[13])
          else 
            # kludge for multiple potentials for same torsion TODO fix!
            puts "WARNING: %9s %8s %8s %8s" %  [tuple, tokens[12], tokens[14], tokens[13]]
            puts "WARNING: %9s %8s %8s %8s" %  [tuple, tors_k[tuple], tors_phase[tuple], tors_n[tuple]]
            tors_k[tuple] = tokens[12]
            tors_phase[tuple] = tokens[14]
            tors_n[tuple] = tokens[13]
          end
        end
      end
    end
  end
end
file.close

atom_total = 0
atom_labels.each do|label,count|
  puts "#%3d %5s %8s %8s" % [count, label, atom_sigma[label], atom_epsilon[label]]
  atom_total += count
end
puts "#no atom types %d: total atoms %s" % [atom_labels.count, atom_total]
puts
puts "epsilon kcal"
atom_labels.each do|label,count|
  puts "%-5s %8s %8.5f" % [equiv_table[label], atom_epsilon[label], Float(atom_sigma[label])*2]
  atom_total += count
end
puts "lennard epsilon geometric 12 6 x13 kcal all 0.5"
puts "12.0"
puts
puts "coul o14 1/6"
puts "X X 7.0"
puts

bond_total = 0
bond_labels.each do|label,count|
  puts "#%3d %5s %8s %8s" % [count, label, bond_k[label], bond_r[label]]
  bond_total += count
end
puts "#no pots %d: total bonds %d" % [bond_labels.count, bond_total]
puts
puts "harm bond kcal"
bond_labels.each do|label,count|
  tokens =  label.split
  puts "%-2s %-2s %10.3f %8s " % [equiv_table[tokens[0]], equiv_table[tokens[1]], Float(bond_k[label])*2, bond_r[label]]
end
puts

three_total = 0
three_labels.each do|label,count|
  puts "#%3d %10s %8s %8s" % [count, label, three_k[label], three_angle[label]]
  three_total += count
end
puts "#no pots %d: total angles %d" % [three_labels.count, three_total]
puts
puts "three bond kcal"
three_labels.each do|label,count|
  tokens = label.split
  puts "%-2s %-2s %-2s %10.3f %8s " % [equiv_table[tokens[1]], equiv_table[tokens[0]], equiv_table[tokens[2]], Float(three_k[label])*2, three_angle[label]]
end
puts

tors_total = 0
tors_labels.each do|label,count|
  puts "#%3d %9s %8s %8s %8s" % [count, label, tors_k[label], tors_phase[label], tors_n[label]]
  tors_total += count
end
puts "#no pots %d: total dihedrals %d" % [tors_labels.count, tors_total]
puts ""
puts "torsion bond kcal"
tors_labels.each do|label,count|
  tokens = label.split
  puts "%-2s %-2s %-2s %-2s %9s %9s %9s" % [equiv_table[tokens[0]], equiv_table[tokens[1]], equiv_table[tokens[2]], equiv_table[tokens[3]], tors_k[label], tors_n[label], tors_phase[label]]
end
puts

improper_total = 0
improper_labels.each do|label,count|
  puts "#%3d %9s %8s %8s %8s" % [count, label, improper_k[label], improper_phase[label], improper_n[label]]
  improper_total += count
end
puts "#no pots %d: total impropers %d" % [improper_labels.count, improper_total]
puts ""
puts "torsion improper kcal"
improper_total = 0
improper_labels.each do|label,count|
  tokens = label.split
  if ((tokens[0] == tokens[1]) && (tokens[1] == tokens[3]))
    multiplicity = 6
  elsif ((tokens[0] == tokens[1]) || (tokens[1] == tokens[3]) || (tokens[0] == tokens[3]))
    multiplicity = 2
  else
    multiplicity = 1
  end
  k = Float(improper_k[label])/multiplicity
  improper_total = improper_total + multiplicity*count
  # swap 2nd and 3rd atom type
  puts "%-2s %-2s %-2s %-2s %10.7f %9s %9s" % [equiv_table[tokens[0]], equiv_table[tokens[2]], equiv_table[tokens[1]], equiv_table[tokens[3]], k, improper_n[label], improper_phase[label]]
end
puts "#total GULP impropers %d" % [improper_total]

