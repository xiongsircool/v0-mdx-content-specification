import type { ReactNode } from "react"
import { cn } from "@/lib/utils"

interface StepsProps {
  children: ReactNode
}

export function Steps({ children }: StepsProps) {
  return (
    <div className="my-6">
      <ol className="relative space-y-6 [counter-reset:step]">{children}</ol>
    </div>
  )
}

interface StepProps {
  children: ReactNode
  className?: string
}

export function Step({ children, className }: StepProps) {
  return (
    <li
      className={cn(
        "relative flex gap-4 pb-6 [counter-increment:step] last:pb-0",
        "before:absolute before:left-0 before:top-0 before:flex before:h-8 before:w-8 before:items-center before:justify-center",
        "before:rounded-full before:bg-primary before:text-sm before:font-semibold before:text-primary-foreground",
        "before:content-[counter(step)]",
        "after:absolute after:left-4 after:top-8 after:h-full after:w-px after:bg-border last:after:hidden",
        className,
      )}
    >
      <div className="flex-1 pl-12">
        <div className="prose prose-sm max-w-none dark:prose-invert">{children}</div>
      </div>
    </li>
  )
}
