  section .data
val dd 513
; 513 = 0b1000000001
; n = 2
; should have eax = rev[0,1,0,...0] = [0,...0,1,0] = 2
;
; 123 = 0b1111011
; n = 6                                                 
; should have eax = rev[6%2,15%2,15%2,0,...,0]          
;                 = rev[0,1,1,0,...,0]                  
;                 = [0,...,0,1,1,0]                     
;                 = 6

printf_format db "val = %u",0x0a,0

  section .text
  global main
  extern printf

main:
  mov   edi, [val]
  mov   edx,edi
  shr   edx,1
  and   edx,0x55555555
  sub   edi,edx
  mov   eax,edi
  shr   edi,0x2
  and   eax,0x33333333
  and   edi,0x33333333
  add   edi,eax
  mov   eax,edi
  shr   eax,0x4
  add   eax,edi
  and   eax,0x0f0f0f0f
  imul  eax,eax,0x01010101
  shr   eax,0x18
  lea   edi, [printf_format]
  mov   esi, eax
  xor   eax, eax
  call printf
  ret 
