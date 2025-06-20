# MIDI -> my sheet with no argument.

# mido powered by: https://mido.readthedocs.io/en/latest/installing.html

from mido import MidiFile
import sys

# oct = int(input())

# def turn_to_note_name(note):
#   return 1    

midi_sheet = MidiFile('reyichang.mid')

res = open('result.txt','w',encoding='utf-8')

sys.stdout = res

for i, track in enumerate(midi_sheet.tracks):
    print('Track {}: {}'.format(i,track.name))
    for msg in track:
        print(msg)