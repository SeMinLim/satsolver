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
    .sdata : ALIGN(0x20000000) {
	*(.sdata)
    } > ram
    .srodata : {
        *(.srodata)
    } > ram
    .rodata : {
	*(.rodata)
    } > ram
    .riscv.attributes : {
	*(.riscv.attributes)
    } > ram
    .bss : {
        *(.bss)
    } > ram
    .sbss : {
	*(.sbss)
    } > ram
}
