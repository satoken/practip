#!/usr/bin/env ruby

while l=gets
  l=l.chomp.split
  if l.size==4
    puts "{"+["\"#{l[1]}\"", 0.0].join(",")+"},"
  elsif l.size==2
    puts "{"+["\"#{l[0]}\"", l[1]].join(", ")+"},"
  end
end
#puts "{"+["NULL", 0.0].join(",")+"},"
