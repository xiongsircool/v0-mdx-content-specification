import { Badge as UIBadge } from "@/components/ui/badge"
import type { ReactNode } from "react"

interface BadgeProps {
  children: ReactNode
  variant?: "default" | "secondary" | "destructive" | "outline"
}

export function Badge({ children, variant = "secondary" }: BadgeProps) {
  return (
    <UIBadge variant={variant} className="mx-1">
      {children}
    </UIBadge>
  )
}
