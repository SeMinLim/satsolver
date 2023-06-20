MEMORY
{
    ram (rwx) : ORIGIN = 0x00000000, LENGTH = 0xffffffff /* 4GB */
}
SECTIONS
{
    . = 0x0;
    .text : {
        *(.text)
    } > ram
    .rodata : ALIGN(0x4000){
	*(.rodata)
    } > ram
    .sdata : {
	*(.sdata)
    } > ram
    /*.srodata : {
	*(.srodata)
    } > ram*/ /* For printing output result*/
    .sbss : {
	*(.sbss)
    } > ram
    .bss : {
        *(.bss)
    } > ram
    .riscv.attributes : {
	*(.riscv.attributes)
    } > ram
}
