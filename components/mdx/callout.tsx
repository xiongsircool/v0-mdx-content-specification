import type { ReactNode } from "react"
import { cn } from "@/lib/utils"
import { AlertCircle, CheckCircle, AlertTriangle, XCircle } from "lucide-react"

interface CalloutProps {
  type?: "info" | "success" | "warning" | "danger"
  children: ReactNode
}

export function Callout({ type = "info", children }: CalloutProps) {
  const config = {
    info: {
      icon: AlertCircle,
      className: "border-blue-200 bg-blue-50 text-blue-900 dark:border-blue-800 dark:bg-blue-950 dark:text-blue-100",
      iconClassName: "text-blue-600 dark:text-blue-400",
    },
    success: {
      icon: CheckCircle,
      className:
        "border-green-200 bg-green-50 text-green-900 dark:border-green-800 dark:bg-green-950 dark:text-green-100",
      iconClassName: "text-green-600 dark:text-green-400",
    },
    warning: {
      icon: AlertTriangle,
      className:
        "border-yellow-200 bg-yellow-50 text-yellow-900 dark:border-yellow-800 dark:bg-yellow-950 dark:text-yellow-100",
      iconClassName: "text-yellow-600 dark:text-yellow-400",
    },
    danger: {
      icon: XCircle,
      className: "border-red-200 bg-red-50 text-red-900 dark:border-red-800 dark:bg-red-950 dark:text-red-100",
      iconClassName: "text-red-600 dark:text-red-400",
    },
  }

  const { icon: Icon, className, iconClassName } = config[type]

  return (
    <div className={cn("my-6 flex gap-3 rounded-lg border p-4", className)}>
      <Icon className={cn("h-5 w-5 flex-shrink-0 mt-0.5", iconClassName)} />
      <div className="flex-1">{children}</div>
    </div>
  )
}
